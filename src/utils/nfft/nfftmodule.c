#include <Python.h>
#include "structmember.h"
#include <numpy/arrayobject.h>
#include <nfft3.h>
#include <math.h>
#include <stdio.h>

#if defined(ENABLE_THREADS)
#include <omp.h>
#include "nfft3util.h"
#endif



PyDoc_STRVAR(nfft__doc__, "nfft(real_space, coordinates)\n\nCalculate nfft from arbitrary dimensional array.\nreal_space should be an array (or any object that can trivially be converted to one.\ncoordinates should be a NxD array where N is the number of points where the Fourier transform should be evaluated and D is the dimensionality of the input array");
static PyObject *nfft(PyObject *self, PyObject *args, PyObject *kwargs)
{
  PyObject *in_obj, *coord_obj;

  
  
  static char *kwlist[] = {"real_space", "coordinates", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OO", kwlist, &in_obj, &coord_obj)) {
    return NULL;
  }
  
  PyObject *coord_array = PyArray_FROM_OTF(coord_obj, NPY_DOUBLE, NPY_IN_ARRAY);
  PyObject *in_array = PyArray_FROM_OTF(in_obj, NPY_COMPLEX128, NPY_IN_ARRAY);
  if (coord_array == NULL || in_array == NULL) {
    Py_XDECREF(coord_array);
    Py_XDECREF(in_array);
    return NULL;
  }

  int ndim = PyArray_NDIM(in_array);
  if (ndim <= 0) {
    PyErr_SetString(PyExc_ValueError, "Input array can't be 0 dimensional\n");
    return NULL;
  }

  if ((PyArray_NDIM(coord_array) != 2 || PyArray_DIM(coord_array, 1) != ndim) && (ndim != 1 || PyArray_NDIM(coord_array) != 1)) {
    PyErr_SetString(PyExc_ValueError, "Coordinates must be given as array of dimensions [NUMBER_OF_POINTS, NUMBER_OF_DIMENSIONS] of [NUMBER_OF_POINTS for 1D transforms.\n");
    Py_XDECREF(coord_array);
    Py_XDECREF(in_array);
    return NULL;
  }
  int number_of_points = (int) PyArray_DIM(coord_array, 0);

  
  nfft_plan my_plan;
  int total_number_of_pixels = 1;
  int dims[ndim];
  int dim;
  for (dim = 0; dim < ndim; ++dim) {
    dims[dim] = (int)PyArray_DIM(in_array, dim);
    total_number_of_pixels *= dims[dim];
  }

  #if defined(ENABLE_THREADS)
  printf("OMP_NUM_THREADS=%s\n",getenv("OMP_NUM_THREADS"));   
  printf("nthreads = %d\n", nfft_get_num_threads());
  fftw_init_threads();
  #endif

  nfft_init(&my_plan, ndim, dims, number_of_points);
  memcpy(my_plan.f_hat, PyArray_DATA(in_array), total_number_of_pixels*sizeof(fftw_complex));
  memcpy(my_plan.x, PyArray_DATA(coord_array), ndim*number_of_points*sizeof(double));
  
  // As of NFFT 3.3, "nfft_flags" has been renamed to "flags"
  if (my_plan.nfft_flags &PRE_PSI) {
    nfft_precompute_one_psi(&my_plan);
  }

  nfft_trafo(&my_plan);

  int out_dim[] = {number_of_points};
  PyObject *out_array = (PyObject *)PyArray_FromDims(1, out_dim, NPY_COMPLEX128);
  memcpy(PyArray_DATA(out_array), my_plan.f, number_of_points*sizeof(fftw_complex));

  // Clean up memory
  
  nfft_finalize(&my_plan);

  #if defined(ENABLE_THREADS)
  fftw_cleanup_threads();
  #endif

  Py_XDECREF(coord_array);
  Py_XDECREF(in_array);
  
  return out_array;
}

static PyMethodDef NfftMethods[] = {
  {"nfft", (PyCFunction)nfft, METH_VARARGS|METH_KEYWORDS, nfft__doc__},
  {NULL, NULL, 0, NULL}
};

// Macro to make module definition compatible with python 2 and 3
// taken from: http://python3porting.com/cextensions.html#module-initialization
#if PY_MAJOR_VERSION >= 3
#define MOD_ERROR_VAL NULL
#define MOD_SUCCESS_VAL(val) val
#define MOD_INIT(name) PyMODINIT_FUNC PyInit_##name(void)
#define MOD_DEF(ob, name, doc, methods) \
  static struct PyModuleDef moduledef = { \
    PyModuleDef_HEAD_INIT, name, doc, -1, methods, }; \
  ob = PyModule_Create(&moduledef);
#else
#define MOD_ERROR_VAL
#define MOD_SUCCESS_VAL(val)
#define MOD_INIT(name) void init##name(void)
#define MOD_DEF(ob, name, doc, methods) \
  ob = Py_InitModule3(name, methods, doc);
#endif

MOD_INIT(nfft)
{
  import_array();
  PyObject *m;
  MOD_DEF(m, "nfft", "Nonequispaced FFT tools.", NfftMethods)
  if (m == NULL)
    return MOD_ERROR_VAL;
  return MOD_SUCCESS_VAL(m);  
}
