#include <Python.h>
#include "structmember.h"
#include <numpy/arrayobject.h>
#include <nfft3.h>
#include <math.h>
#include <stdio.h>

//#include "nfftclassmodule.h"

static void nfft_1d_carrays(double *in, double *coord, double *out, int number_of_pixels, int number_of_points)
{
  nfft_plan my_plan;
  nfft_init_1d(&my_plan, number_of_pixels, number_of_points);

  memcpy(my_plan.x, coord, number_of_points*sizeof(double));
  //my_plan.x = coord
  memcpy(my_plan.f_hat, in, 2*number_of_pixels*sizeof(double));

  if (my_plan.nfft_flags &PRE_PSI) {
    nfft_precompute_one_psi(&my_plan);
  }
  
  nfft_trafo(&my_plan);

  memcpy(out, my_plan.f, 2*number_of_points*sizeof(double));

  nfft_finalize(&my_plan);
}

static void nfft_3d_carrays(double *in, double *coord, double *out, int number_of_pixels_x,
			    int number_of_pixels_y, int number_of_pixels_z, int number_of_points)
{
  nfft_plan my_plan;
  nfft_init_3d(&my_plan, number_of_pixels_z, number_of_pixels_y, number_of_pixels_x, number_of_points);

  memcpy(my_plan.x, coord, 3*number_of_points*sizeof(double));
  memcpy(my_plan.f_hat, in, 2*number_of_pixels_x*number_of_pixels_y*number_of_pixels_z*sizeof(double));

  if (my_plan.nfft_flags &PRE_PSI) {
    nfft_precompute_one_psi(&my_plan);
  }
  
  nfft_trafo(&my_plan);
  memcpy(out, my_plan.f, 2*number_of_points*sizeof(double));

  nfft_finalize(&my_plan);
}

static PyObject *nfft(PyObject *self, PyObject *args)
{
  PyObject *vecin_obj, *veccoord_obj;
  double *cin, *cout, *ccoord;

  int number_of_pixels, number_of_points;

  if (!PyArg_ParseTuple(args, "OO", &vecin_obj, &veccoord_obj)) {
    return NULL;
  }

  PyObject *veccoord_array = PyArray_FROM_OTF(veccoord_obj, NPY_DOUBLE, NPY_IN_ARRAY);
  PyObject *vecin_array = PyArray_FROM_OTF(vecin_obj, NPY_COMPLEX128, NPY_IN_ARRAY);
  if (veccoord_array == NULL || vecin_array == NULL) {
    Py_XDECREF(veccoord_array);
    Py_XDECREF(vecin_array);
    return NULL;
  }

  if (PyArray_NDIM(veccoord_array) != 1 || PyArray_NDIM(vecin_array) != 1) {
    PyErr_SetString(PyExc_ValueError, "Input must all be 1D arrays.");
    Py_XDECREF(veccoord_array);
    Py_XDECREF(vecin_array);
    return NULL;
  }
  number_of_pixels = (int)PyArray_DIM(vecin_array, 0);
  number_of_points = (int)PyArray_DIM(veccoord_array, 0);

  int vecout_dims[] = {number_of_points};
  PyObject *vecout_array = PyArray_FromDims(1, vecout_dims, NPY_COMPLEX128);

  ccoord = (double *)PyArray_DATA(veccoord_array);
  cin = (double *)PyArray_DATA(vecin_array);
  cout = (double *)PyArray_DATA(vecout_array);

  nfft_1d_carrays(cin, ccoord, cout, number_of_pixels, number_of_points);
  
  Py_XDECREF(veccoord_array);
  Py_XDECREF(vecin_array);
  return vecout_array;
}

static PyObject *nfft_inplace(PyObject *self, PyObject *args)
{
  PyObject *in_obj, *out_obj, *coord_obj;
  double *cin, *cout, *ccoord;

  int number_of_pixels, number_of_points;

  if (!PyArg_ParseTuple(args, "OOO", &in_obj, &coord_obj, &out_obj)) {
    return NULL;
  }

  PyObject *coord_array = PyArray_FROM_OTF(coord_obj, NPY_DOUBLE, NPY_IN_ARRAY);
  PyObject *in_array = PyArray_FROM_OTF(in_obj, NPY_COMPLEX128, NPY_IN_ARRAY);

  if (!PyArray_Check(out_obj)) {
    PyErr_SetString(PyExc_ValueError, "Output must be numpy.array of dtype complex128");
    return NULL;
  }

  if (PyArray_DESCR(out_obj)->type_num != NPY_COMPLEX128) {
    PyErr_SetString(PyExc_ValueError, "Output array must be of type complex128");
    Py_XDECREF(coord_array);
    Py_XDECREF(in_array);
    return NULL;
  }
  PyObject *out_array = PyArray_FROM_OTF(out_obj, NPY_COMPLEX128, NPY_INOUT_ARRAY);

  if (coord_array == NULL || in_array == NULL || out_array == NULL)
    goto fail;
  
  if (PyArray_NDIM(coord_array) != 1 || PyArray_NDIM(in_array) != 1 || PyArray_NDIM(out_array) != 1) {
    PyErr_SetString(PyExc_ValueError, "Input must all be 1D arrays.");
    goto fail;
  }
  number_of_pixels = (int)PyArray_DIM(in_array, 0);
  number_of_points = (int)PyArray_DIM(out_array, 0);
  ccoord = (double *)PyArray_DATA(coord_array);
  cin = (double *)PyArray_DATA(in_array);
  cout = (double *)PyArray_DATA(out_array);
  if ((int)PyArray_DIM(coord_array, 0) != number_of_points) {
    PyErr_SetString(PyExc_ValueError, "Coordinates and output must be the same length.");
    goto fail;
  }

  nfft_1d_carrays(cin, ccoord, cout, number_of_pixels, number_of_points);

  Py_XDECREF(coord_array);
  Py_XDECREF(in_array);
  Py_XDECREF(out_array);
  return Py_BuildValue("i", 1);

 fail:
  Py_XDECREF(coord_array);
  Py_XDECREF(in_array);
  Py_XDECREF(out_array);
  return NULL;
}

static PyObject *nfft3(PyObject *self, PyObject *args)
{
  PyObject *in_obj, *coord_obj;
  double *cin, *cout, *ccoord;

  if (!PyArg_ParseTuple(args, "OO", &in_obj, &coord_obj)) {
    return NULL;
  }

  PyObject *coord_array = PyArray_FROM_OTF(coord_obj, NPY_DOUBLE, NPY_IN_ARRAY);
  PyObject *in_array = PyArray_FROM_OTF(in_obj, NPY_COMPLEX128, NPY_IN_ARRAY);
  if (coord_array == NULL || in_array == NULL) {
    Py_XDECREF(coord_array);
    Py_XDECREF(in_array);
    return NULL;
  }

  if (PyArray_NDIM(coord_array) != 2 || PyArray_DIM(coord_array, 1) != 3) {
    PyErr_SetString(PyExc_ValueError, "Coordinate array must be size Nx3");
    Py_XDECREF(coord_array);
    Py_XDECREF(in_array);
    return NULL;
  }

  if (PyArray_NDIM(in_array) != 3) {
    PyErr_SetString(PyExc_ValueError, "Input must be a 3D array.");
    Py_XDECREF(coord_array);
    Py_XDECREF(in_array);
    return NULL;
  }

  int number_of_pixels_z = (int)PyArray_DIM(in_array, 0);
  int number_of_pixels_y = (int)PyArray_DIM(in_array, 1);
  int number_of_pixels_x = (int)PyArray_DIM(in_array, 2);
  int number_of_points = (int)PyArray_DIM(coord_array, 0);

  int out_dims[] = {number_of_points};
  PyObject *out_array = PyArray_FromDims(1, out_dims, NPY_COMPLEX128);

  ccoord = (double *)PyArray_DATA(coord_array);
  cin = (double *)PyArray_DATA(in_array);
  cout = (double *)PyArray_DATA(out_array);

  nfft_3d_carrays(cin, ccoord, cout, number_of_pixels_x, number_of_pixels_y, number_of_pixels_z, number_of_points);

  Py_XDECREF(coord_array);
  Py_XDECREF(in_array);
  return out_array;
}

static PyObject *nfft3_inplace(PyObject *self, PyObject *args)
{
  PyObject *in_obj, *coord_obj, *out_obj;
  double *cin, *cout, *ccoord;

  if (!PyArg_ParseTuple(args, "OOO", &in_obj, &coord_obj, &out_obj)) {
    return NULL;
  }

  PyObject *coord_array = PyArray_FROM_OTF(coord_obj, NPY_DOUBLE, NPY_IN_ARRAY);
  PyObject *in_array = PyArray_FROM_OTF(in_obj, NPY_COMPLEX128, NPY_IN_ARRAY);
  if (coord_array == NULL || in_array == NULL) {
    Py_XDECREF(coord_array);
    Py_XDECREF(in_array);
    return NULL;
  }

  if (PyArray_NDIM(coord_array) != 2 || PyArray_DIM(coord_array, 1) != 3) {
    PyErr_SetString(PyExc_ValueError, "Coordinate array must be size Nx3");
    Py_XDECREF(coord_array);
    Py_XDECREF(in_array);
    return NULL;
  }

  if (PyArray_NDIM(in_array) != 3) {
    PyErr_SetString(PyExc_ValueError, "Input must be a 3D array.");
    Py_XDECREF(coord_array);
    Py_XDECREF(in_array);
    return NULL;
  }

  int number_of_pixels_z = (int)PyArray_DIM(in_array, 0);
  int number_of_pixels_y = (int)PyArray_DIM(in_array, 1);
  int number_of_pixels_x = (int)PyArray_DIM(in_array, 2);
  int number_of_points = (int)PyArray_DIM(coord_array, 0);

  if (!PyArray_Check(out_obj)) {
    PyErr_SetString(PyExc_ValueError, "Output must be numpy.array of dtype complex128");
    return NULL;
  }

  if (PyArray_DESCR(out_obj)->type_num != NPY_COMPLEX128) {
    PyErr_SetString(PyExc_ValueError, "Output array must be of dtype complex128");
    Py_XDECREF(coord_array);
    Py_XDECREF(in_array);
    return NULL;
  }
  PyObject *out_array = PyArray_FROM_OTF(out_obj, NPY_COMPLEX128, NPY_INOUT_ARRAY);

  if (PyArray_NDIM(out_array) != 1 || PyArray_DIM(out_array, 0) != number_of_points) {
    PyErr_SetString(PyExc_ValueError, "Output array must be one dimensional and same length as coordinates array.");
    Py_XDECREF(coord_array);
    Py_XDECREF(in_array);
    Py_XDECREF(out_array);
    return NULL;
  }

  ccoord = (double *)PyArray_DATA(coord_array);
  cin = (double *)PyArray_DATA(in_array);
  cout = (double *)PyArray_DATA(out_array);

  nfft_3d_carrays(cin, ccoord, cout, number_of_pixels_x, number_of_pixels_y, number_of_pixels_z, number_of_points);

  Py_XDECREF(coord_array);
  Py_XDECREF(in_array);
  Py_XDECREF(out_array);
  return Py_BuildValue("i", 1);
}


static PyObject *nfftn(PyObject *self, PyObject *args)
{
  PyObject *in_obj, *coord_obj;

  if (!PyArg_ParseTuple(args, "OO", &in_obj, &coord_obj)) {
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
  int M_total = 1;
  int dims[ndim];
  for (int dim = 0; dim < ndim; ++dim) {
    dims[dim] = (int)PyArray_DIM(in_array, dim);
    M_total *= dims[dim];
  }
  nfft_init(&my_plan, ndim, dims, M_total);
  memcpy(my_plan.f_hat, PyArray_DATA(in_array), M_total*sizeof(fftw_complex));
  memcpy(my_plan.x, PyArray_DATA(coord_array), ndim*number_of_points*sizeof(double));
  
  if (my_plan.nfft_flags &PRE_PSI) {
    nfft_precompute_one_psi(&my_plan);
  }

  nfft_trafo(&my_plan);
  
  //npy_intp out_dim[] = {number_of_points};
  int out_dim[] = {number_of_points};
  //PyObject *vecout_array = PyArray_FromDims(1, vecout_dims, NPY_COMPLEX128);
  //PyObject *out_array = (PyObject *)PyArray_SimpleNew(1, out_dim, NPY_COMPLEX128);
  PyObject *out_array = (PyObject *)PyArray_FromDims(1, out_dim, NPY_COMPLEX128);
  memcpy(PyArray_DATA(out_array), my_plan.f, number_of_points*sizeof(fftw_complex));

  nfft_finalize(&my_plan);
  return out_array;
}

static PyObject *nfftn_inplace(PyObject *self, PyObject *args)
{
  PyObject *in_obj, *coord_obj, *out_obj;

  if (!PyArg_ParseTuple(args, "OOO", &in_obj, &coord_obj, &out_obj)) {
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

  if (!PyArray_Check(out_obj)) {
    PyErr_SetString(PyExc_ValueError, "Output must be numpy.array of dtype complex128");
    Py_XDECREF(coord_array);
    Py_XDECREF(in_array);
    return NULL;
  }

  if (PyArray_DESCR(out_obj)->type_num != NPY_COMPLEX128) {
    PyErr_SetString(PyExc_ValueError, "Output array must be of dtype complex128");
    Py_XDECREF(coord_array);
    Py_XDECREF(in_array);
    return NULL;
  }

  PyObject *out_array = PyArray_FROM_OTF(out_obj, NPY_COMPLEX128, NPY_INOUT_ARRAY);
  if (PyArray_NDIM(out_array) != 1 || PyArray_DIM(out_array, 0) != number_of_points) {
    PyErr_SetString(PyExc_ValueError,  "Output array must be one dimensional and same length as coordinates array.");
    Py_XDECREF(coord_array);
    Py_XDECREF(in_array);
    Py_XDECREF(out_array);
    return NULL;
  }

  nfft_plan my_plan;
  int M_total = 1;
  int dims[ndim];
  for (int dim = 0; dim < ndim; ++dim) {
    dims[dim] = (int)PyArray_DIM(in_array, dim);
    M_total *= dims[dim];
  }
  nfft_init(&my_plan, ndim, dims, M_total);
  memcpy(my_plan.f_hat, PyArray_DATA(in_array), M_total*sizeof(fftw_complex));
  memcpy(my_plan.x, PyArray_DATA(coord_array), ndim*number_of_points*sizeof(double));
  
  if (my_plan.nfft_flags &PRE_PSI) {
    nfft_precompute_one_psi(&my_plan);
  }

  nfft_trafo(&my_plan);
  
  memcpy(PyArray_DATA(out_array), my_plan.f, number_of_points*sizeof(fftw_complex));

  nfft_finalize(&my_plan);
  return Py_BuildValue("i", 1);
}


/*
static PyObject *nfft3(PyObject *self, PyObject *args)
{
  PyObject *vecin_obj, *veccoord_obj;
  double *cin, *cout, *ccoord;

  int number_of_pixels_x, number_of_pixels_y, number_of_pixels_z;
  int number_of_points;
  int i;

  if (!PyArg_ParseTuple(args, "OO", &vecin_obj, &veccoord_obj))
    return NULL;
  PyObject *veccoord_array = PyArray_FROM_OTF(veccoord_obj, NPY_DOUBLE, NPY_IN_ARRAY);
  PyObject *vecin_array = PyArray_FROM_OTF(vecin_obj, NPY_COMPLEXLTR, NPY_IN_ARRAY);
  if (veccoord_array == NULL || vecin_array == NULL) {
    Py_XDECREF(veccoord_array);
    Py_XDECREF(vecin_array);
    return NULL;
  }

  if (PyArray_NDIM(vecin_array) != 3) {
    Py_XDECREF(veccoord_array);
    Py_XDECREF(vecin_array);
    PyErr_SetString(PyExc_ValueError, "Input array must be three dimensional");
    return NULL;
  }
  number_of_pixels_z = (int)PyArray_DIM(vecin_array, 0);
  number_of_pixels_y = (int)PyArray_DIM(vecin_array, 1);
  number_of_pixels_x = (int)PyArray_DIM(vecin_array, 2);

  if (PyArray_NDIM(veccoord_array) != 2 || PyArray_DIM(veccoord_array, 1) != 3) {
    Py_XDECREF(veccoord_array);
    Py_XDECREF(vecin_array);
    PyErr_SetString(PyExc_ValueError, "Coordinate array must be three dimensional");
    return NULL;
  }    
  number_of_points = (int)PyArray_DIM(veccoord_array, 0);

  ccoord = PyArray_DATA(veccoord_array);
  cin = PyArray_DATA(vecin_array);

  int vecout_dims[] = {number_of_points};
  PyObject *vecout_array = PyArray_FromDims(1, vecout_dims, NPY_COMPLEX128);
  cout = PyArray_DATA(vecout_array);



  //return Py_BuildValue("i", 1);
  Py_XDECREF(veccoord_array);
  Py_XDECREF(vecin_array);
  return vecout_array;
}
*/

typedef struct {
  PyObject_HEAD
  nfft_plan my_plan;
  PyArrayObject *real_map;
  fftw_complex *sneaky_ref;
  int ndim;
}Transformer;

static int Transformer_init(Transformer *self, PyObject *args, PyObject *kwds)
{
  PyObject *input_obj;
  if (!PyArg_ParseTuple(args, "O", &input_obj)) {
    return -1;
  }
  PyObject *input_array = PyArray_FROM_OTF(input_obj, NPY_COMPLEX128, NPY_IN_ARRAY);
  if (input_array == NULL) {
    return -1;
  }
  self->real_map = (PyArrayObject *)input_array;

  self->ndim = PyArray_NDIM(input_array);
  int M_total = 1;
  int dims[self->ndim];
  for (int dim = 0; dim < self->ndim; ++dim) {
    dims[dim] = (int) PyArray_DIM(input_array, dim);
    M_total *= dims[dim];
  }
  nfft_init(&self->my_plan, self->ndim, dims, M_total);
  //free(self->my_plan.f_hat);
  self->sneaky_ref = self->my_plan.f_hat;
  self->my_plan.f_hat = (fftw_complex *)PyArray_DATA(self->real_map);
  return 0;
}

static void Transformer_dealloc(Transformer *self)
{
  Py_XDECREF(self->real_map);
  self->my_plan.f_hat = self->sneaky_ref;
  nfft_finalize(&self->my_plan);
  self->ob_type->tp_free((PyObject *)self);
}

static PyObject *Transformer_transform(Transformer *self, PyObject *args)
{
  PyObject *input_obj;
  if (!PyArg_ParseTuple(args, "O", &input_obj)) {
    return NULL;
  }
  PyObject *coordinates_array = PyArray_FROM_OTF(input_obj, NPY_FLOAT64, NPY_IN_ARRAY);

  if (coordinates_array == NULL) {
    return NULL;
  }
  if ((PyArray_NDIM(coordinates_array) != 2 || PyArray_DIM(coordinates_array, 1) != self->ndim) && (self->ndim != 1 || PyArray_NDIM(coordinates_array) != 1)) {
    PyErr_SetString(PyExc_ValueError, "Coordinates must be given as array of dimensions [NUMBER_OF_POINTS, NUMBER_OF_DIMENSIONS] of [NUMBER_OF_POINTS for 1D transforms.\n");
    Py_XDECREF(coordinates_array);
    return NULL;
  }
  int number_of_points = (int) PyArray_DIM(coordinates_array, 0);
  
  memcpy(self->my_plan.x, PyArray_DATA(coordinates_array), self->ndim*number_of_points*sizeof(double));
  
  if (self->my_plan.nfft_flags &PRE_PSI) {
    nfft_precompute_one_psi(&self->my_plan);
  }

  nfft_trafo(&self->my_plan);
  
  npy_intp out_dim[] = {number_of_points};
  PyObject *out_array = (PyObject *)PyArray_SimpleNew(1, out_dim, NPY_COMPLEX128);
  memcpy(PyArray_DATA(out_array), self->my_plan.f, number_of_points*sizeof(fftw_complex));
  return out_array;
  /*
  nfft_plan my_plan;
  nfft_init_3d(&my_plan, number_of_pixels_z, number_of_pixels_y, number_of_pixels_x, number_of_points);

  memcpy(my_plan.x, coord, 3*number_of_points*sizeof(double));
  memcpy(my_plan.f_hat, in, 2*number_of_pixels_x*number_of_pixels_y*number_of_pixels_z*sizeof(double));

  if (my_plan.nfft_flags &PRE_PSI) {
    nfft_precompute_one_psi(&my_plan);
  }
  
  nfft_trafo(&my_plan);
  memcpy(out, my_plan.f, 2*number_of_points*sizeof(double));

  nfft_finalize(&my_plan);
  */

  return Py_BuildValue("");
}

static PyObject *ndim(Transformer *self, PyObject *args, PyObject *kwds)
{
  if (self->ndim > 0) {
    return Py_BuildValue("i", self->ndim);
  } else {
    return Py_BuildValue("");
  }
}

static PyMemberDef Transformer_members[] = {
  {"real_map", T_OBJECT_EX, offsetof(Transformer, real_map), 0, "Real space map."},
  {NULL}
};

static PyMethodDef Transformer_methods[] = {
  {"ndim", (PyCFunction) ndim, METH_VARARGS, "ndim()\n\nGet the number of dimensions."},
  {"transform", (PyCFunction) Transformer_transform, METH_VARARGS, "transform(coordinates)\n\nReturns transformation at given coordintase.\n\nCoordinate array has format [NUMBER_OF_POINTS, NUMBER_OF_DIMENSIONS] and can be of any type that can direcly be converted to an ndarray."},
  {NULL}
};

static PyTypeObject TransformerType = {
   PyObject_HEAD_INIT(NULL)
   0,                         /* ob_size */
   "Transformer",         /* tp_name */
   sizeof(Transformer),   /* tp_basicsize */
   0,                         /* tp_itemsize */
   (destructor)Transformer_dealloc, /* tp_dealloc */
   0,                         /* tp_print */
   0,                         /* tp_getattr */
   0,                         /* tp_setattr */
   0,                         /* tp_compare */
   0,                         /* tp_repr */
   0,                         /* tp_as_number */
   0,                         /* tp_as_sequence */
   0,                         /* tp_as_mapping */
   0,                         /* tp_hash */
   0,                         /* tp_call */
   0,                         /* tp_str */
   0,                         /* tp_getattro */
   0,                         /* tp_setattro */
   0,                         /* tp_as_buffer */
   Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /* tp_flags*/
   "Transformer(real_space)\n\nCreates an object that can be used to calculate multiple nfft transforms from the same array.",  /* tp_doc */
   0,                         /* tp_traverse */
   0,                         /* tp_clear */
   0,                         /* tp_richcompare */
   0,                         /* tp_weaklistoffset */
   0,                         /* tp_iter */
   0,                         /* tp_iternext */
   Transformer_methods,   /* tp_methods */
   Transformer_members,   /* tp_members */
   0,                         /* tp_getset */
   0,                         /* tp_base */
   0,                         /* tp_dict */
   0,                         /* tp_descr_get */
   0,                         /* tp_descr_set */
   0,                         /* tp_dictoffset */
   (initproc)Transformer_init,  /* tp_init */
   0,                         /* tp_alloc */
   0,                         /* tp_new */
};

static PyMethodDef NfftMethods[] = {
  {"nfft", nfft, METH_VARARGS, "nfft(real_space, coordinates)\n\nCalculate 1d nfft.\n\real_space should be 1d array.\ncoordinates should be a 1d array of the coordinates where the Fourier transform should be evaluated."},
  {"nfft_inplace", nfft_inplace, METH_VARARGS, "nfft_inplace(real_space, coordinates, output_array)\n\nCalculate 1d nfft.\n\real space should be 1d array.\ncoordinates should be a 1d array of the coordinates where the Fourier transform should be evaluated\noutput_array should be ndarray of type complex128. The is written to here, if the array is a continuous block in memory this can speed up the calculation."},
  {"nfft3", nfft3, METH_VARARGS, "nfft3(real_space, coordinates)\n\nCalculate 3d nfft.\n\real_space should be 3d array.\ncoordinates should be a Nx3 array where N is the number of points where the Fourier transform should be evaluated."},
  {"nfft3_inplace", nfft3_inplace, METH_VARARGS, "nfft3(real_space, coordinates)\n\nCalculate 3d nfft.\n\real_space should be 3d array.\ncoordinates should be a Nx3 array where N is the number of points where the Fourier transform should be evaluated\noutput_array should be ndarray of type complex128. The is written to here, if the array is a continuous block in memory this can speed up the calculation."},
  {"nfftn", nfftn, METH_VARARGS, "nfft3(real_space, coordinates)\n\nCalculate nfft from arbitrary dimensional array.\n\real_space should be an array (or any object that can trivially be converted to one.\ncoordinates should be a NxD array where N is the number of points where the Fourier transform should be evaluated and D is the dimensionality of the input array"},
  {"nfftn_inplace", nfftn_inplace , METH_VARARGS, "nfft3(real_space, coordinates)\n\nCalculate nfft from arbitrary dimensional array.\n\real_space should be an array (or any object that can trivially be converted to one.\ncoordinates should be a NxD array where N is the number of points where the Fourier transform should be evaluated and D is the dimensionality of the input array\noutput_array should be ndarray of type complex128. The is written to here, if the array is a continuous block in memory this can speed up the calculation."},
  {NULL, NULL, 0, NULL}
};



PyMODINIT_FUNC initnfft(void)
{
  import_array();
  PyObject *m = Py_InitModule3("nfft", NfftMethods, "Nonequispaced FFT tools.");
  if (m == NULL)
    return;

  TransformerType.tp_new = PyType_GenericNew;
  if (PyType_Ready(&TransformerType) < 0)
    return;

  Py_INCREF(&TransformerType);
  PyModule_AddObject(m, "Transformer", (PyObject *)&TransformerType);
    
}
/*
int main(int argc, char *argv[])
{
  Py_SetProgramName(argv[0]);
  Py_Initialize();
  initnfft();
  return 0;
}
*/
