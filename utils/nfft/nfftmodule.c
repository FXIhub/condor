#include <Python.h>
#include <numpy/arrayobject.h>
#include <nfft3.h>
#include <math.h>

//static PyObject *NfftError;

PyArrayObject *pyvector(PyObject *objin);
double *pyvector_to_Carrayptrs(PyArrayObject *arrayin);
int not_doublevector(PyArrayObject *vec);

PyArrayObject *pyvector(PyObject *objin) {
  return (PyArrayObject *) PyArray_ContiguousFromObject(objin, NPY_DOUBLE, 1, 1);
}

double *pyvector_to_Carrayptrs(PyArrayObject *arrayin) {
  int n;

  n = arrayin->dimensions[0];
  return (double *) arrayin->data;
}

double *pyvector_to_Carrayptrs_2d(PyArrayObject *arrayin){ //, int *size_x, int *size_y) {
  if (arrayin->nd != 2) {
    PyErr_SetString(PyExc_ValueError, "In pyvector_to_Carrayptrs_2d: array must 2 dimensional.");
    return NULL;
  }

  //size_x = dimensions[0]; size_y = dimensions[1];

  return (double *) arrayin->data;
}

double *pyvector_to_Carrayptrs_3d(PyArrayObject *arrayin) {//, int *size_x, int *size_y, int *size_z) {
  if (arrayin->nd != 3) {
    PyErr_SetString(PyExc_ValueError, "In pyvector_to_Carrayptrs_3d: array must 3 dimensional.");
    return NULL;
  }

  //size_x = dimensions[0]; size_y = dimensions[1]; size_z = dimensions[2];
  return (double *) arrayin->data;
}

int not_doublevector(PyArrayObject *vec) {
  if (vec->descr->type_num != NPY_DOUBLE) {
    PyErr_SetString(PyExc_ValueError, "In not_doublevector: array must be of type Float and 1 dimensional (n).");
    return 1;
  }
  return 0;
}

static PyObject *nfft(PyObject *M, PyObject *args)
{
  PyArrayObject *vecin, *vecout, *veccoord;
  double *cin, *cout, *ccoord;

  int number_of_pixels, number_of_points;
  int i;

  if (!PyArg_ParseTuple(args, "O!O!O!", &PyArray_Type, &veccoord, &PyArray_Type, &vecin, &PyArray_Type, &vecout))
    return NULL;
  if (NULL == veccoord) return NULL;
  if (NULL == vecin) return NULL;
  if (NULL == vecout) return NULL;

  if (not_doublevector(veccoord)) return NULL;
  if (not_doublevector(vecin)) return NULL;
  if (not_doublevector(vecout)) return NULL;

  ccoord = pyvector_to_Carrayptrs(veccoord);
  cin = pyvector_to_Carrayptrs(vecin);
  cout = pyvector_to_Carrayptrs(vecout);

  number_of_pixels = vecin->dimensions[0]/2;
  // printf("number_of_pixels = %d\n", number_of_pixels);
  number_of_points = vecout->dimensions[0]/2;
  // printf("number_of_points = %d\n", number_of_points);
  // printf("veccoord->dimensions[0] = %d\n", (int)veccoord->dimensions[0]);
  if (veccoord->dimensions[0] != number_of_points) return NULL;

  nfft_plan my_plan;
  nfft_init_1d(&my_plan, number_of_pixels, number_of_points);

  for (i = 0; i < number_of_points; ++i) {
    my_plan.x[i] = ccoord[i];
  }

  for (i = 0; i < number_of_pixels; ++i) {
    my_plan.f_hat[i][0] = cin[2*i];
    my_plan.f_hat[i][1] = cin[2*i+1];
  }

  if (my_plan.nfft_flags &PRE_PSI) {
    nfft_precompute_one_psi(&my_plan);
  }
  
  nfft_trafo(&my_plan);

  for (i = 0; i < number_of_points; ++i) {
    cout[2*i] = my_plan.f[i][0];
    cout[2*i+1] = my_plan.f[i][1];
  }

  nfft_finalize(&my_plan);
  
  return Py_BuildValue("i", 1);
}

static PyObject *nfft3(PyObject *self, PyObject *args)
{
  PyArrayObject *vecin_real,*vecin_imag, *vecout_real,*vecout_imag, *veccoord;
  double *cin_real, *cout_real,*cin_imag, *cout_imag, *ccoord;

  int number_of_pixels_x, number_of_pixels_y, number_of_pixels_z;
  int number_of_points;
  int i;

  //printf("init is done\n");

  if (!PyArg_ParseTuple(args, "O!O!O!O!O!", &PyArray_Type, &veccoord, &PyArray_Type, &vecin_real,&PyArray_Type, &vecin_imag, &PyArray_Type, &vecout_real,&PyArray_Type, &vecout_imag))
    return NULL;
  // printf("parsing is done\n");
  if (NULL == veccoord) return NULL;
  if (NULL == vecin_real) return NULL;
  if (NULL == vecin_imag) return NULL;
  if (NULL == vecout_real) return NULL;
  if (NULL == vecout_imag) return NULL;

  // printf("null checking is done\n");

  if (not_doublevector(veccoord)) return NULL;
  if (not_doublevector(vecin_real)) return NULL;
  if (not_doublevector(vecin_imag)) return NULL;
  if (not_doublevector(vecout_real)) return NULL;
  if (not_doublevector(vecout_imag)) return NULL;

  // printf("reading is done\n");


  // printf("converting is done\n");

  number_of_pixels_x = vecin_real->dimensions[0];
  number_of_pixels_y = vecin_real->dimensions[1];
  number_of_pixels_z = vecin_real->dimensions[2];
  //printf("number_of_pixels_x = %d\n", number_of_pixels_x);
  //printf("number_of_pixels_y = %d\n", number_of_pixels_y);
  //printf("number_of_pixels_z = %d\n", number_of_pixels_z);
  //printf("veccoord->dimensions[0] = %d\n", (int)veccoord->dimensions[0]);
  number_of_points = veccoord->dimensions[0]/3;

  //number_of_rpixels_z = veccoord->dimensions[2];
  
  //printf("number_of_points = %d\n", number_of_points);
  //printf("vecout->dimensions[0] = %d\n", (int)vecout_real->dimensions[0]);
  if (veccoord->dimensions[0] != number_of_points*3) return NULL;

  //printf("dimensions checked\n");

  ccoord = pyvector_to_Carrayptrs(veccoord);
  cin_real = pyvector_to_Carrayptrs_3d(vecin_real);
  cin_imag = pyvector_to_Carrayptrs_3d(vecin_imag);
  cout_real = pyvector_to_Carrayptrs(vecout_real);
  cout_imag = pyvector_to_Carrayptrs(vecout_imag);

  nfft_plan my_plan;
  nfft_init_3d(&my_plan, number_of_pixels_x, number_of_pixels_y, number_of_pixels_z, number_of_points);
  //printf("plan initialized\n");

  for (i = 0; i < number_of_points*3; ++i) {
    my_plan.x[i] = ccoord[i];
  }
  //printf("x set\n");

  for (i = 0; i < number_of_pixels_x*number_of_pixels_y*number_of_pixels_z; ++i) {
    my_plan.f_hat[i][0] = cin_real[i];
    my_plan.f_hat[i][1] = cin_imag[i];
  }
  //printf("f_hat set\n");

  if (my_plan.nfft_flags &PRE_PSI) {
    nfft_precompute_one_psi(&my_plan);
    //printf("precomputation done\n");
  }
  //printf("precomputation passed\n");
  
  nfft_trafo(&my_plan);
  //printf("transform done\n");

  for (i = 0; i < number_of_points; ++i) {
    cout_real[i] = my_plan.f[i][0];
    cout_imag[i] = my_plan.f[i][1];
  }
  //printf("result copied\n");

  nfft_finalize(&my_plan);

  return Py_BuildValue("i", 1);
}

static PyMethodDef NfftMethods[] = {
  {"nfft", nfft, METH_VARARGS, "1d nfft."},
  {"nfft3", nfft3, METH_VARARGS, "3d nfft."},
  {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC initnfft_c(void)
{
  (void) Py_InitModule("nfft_c", NfftMethods);
  import_array()
}

int main(int argc, char *argv[])
{
  Py_SetProgramName(argv[0]);
  Py_Initialize();
  initnfft_c();
  return 0;
}
