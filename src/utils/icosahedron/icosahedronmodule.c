#include <Python.h>
#include <numpy/arrayobject.h>
#include <math.h>
#include <stdio.h>

PyDoc_STRVAR(icosahedron__doc__, "icosahedron(array_side, radius, rotation)\n\nGenerate an icosahedron. Radius is given in pixels and is defined as the distance to the corners. Rotation should be tuple of quaternions (w,x,y,z)");
static PyObject *icosahedron(PyObject *self, PyObject *args, PyObject *kwargs)
{
  int image_side = 0.;
  double radius = 0.;
  PyObject *rotation_obj = NULL;
  PyObject *rotation_sequence;

  static char *kwlist[] = {"array_side", "radius", "rotation", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "id|O", kwlist, &image_side, &radius, &rotation_obj)) {
    return NULL;
  }

  if (image_side <= 0) {
    PyErr_SetString(PyExc_ValueError, "Image side must be > 0.");
    return NULL;
  }

  if (radius <= 0.) {
    PyErr_SetString(PyExc_ValueError, "Radius must be > 0.");
  }

  double rotm_11, rotm_12, rotm_13;
  double rotm_21, rotm_22, rotm_23;
  double rotm_31, rotm_32, rotm_33;
  int rotate;
  if (rotation_obj == NULL) {
    rotate = 0;
  } else {
    rotate = 1;
    rotation_sequence = PySequence_Fast(rotation_obj, "Expected a sequence");
    
    long length = PySequence_Size(rotation_sequence);
    if (length != 3 && length != 4) {
      PyErr_SetString(PyExc_ValueError, "Rotation must be of length 4 (quaternion)");
      return NULL;
    }

    PyObject *seq_1 = PySequence_Fast_GET_ITEM(rotation_sequence, 0);
    PyObject *seq_2 = PySequence_Fast_GET_ITEM(rotation_sequence, 1);
    PyObject *seq_3 = PySequence_Fast_GET_ITEM(rotation_sequence, 2);
    PyObject *seq_4 = PySequence_Fast_GET_ITEM(rotation_sequence, 3);
    double quat_1 = PyFloat_AsDouble(seq_1);
    double quat_2 = PyFloat_AsDouble(seq_2);
    double quat_3 = PyFloat_AsDouble(seq_3);
    double quat_4 = PyFloat_AsDouble(seq_4);
      
    double quaternion_norm = sqrt(pow(quat_1, 2) + pow(quat_2, 2) + pow(quat_3, 2) + pow(quat_4, 2));
    quat_1 /= quaternion_norm;
    quat_2 /= quaternion_norm;
    quat_3 /= quaternion_norm;
    quat_4 /= quaternion_norm;
    
    rotm_11 = quat_1*quat_1 + quat_2*quat_2 - quat_3*quat_3 - quat_4*quat_4;
    rotm_12 = 2.*quat_2*quat_3 - 2.*quat_1*quat_4;
    rotm_13 = 2.*quat_2*quat_4 + 2.*quat_1*quat_3;
      
    rotm_21 = 2.*quat_2*quat_3 + 2.*quat_1*quat_4;
    rotm_22 = quat_1*quat_1 - quat_2*quat_2 + quat_3*quat_3 - quat_4*quat_4;
    rotm_23 = 2.*quat_3*quat_4 - 2.*quat_1*quat_2;

    rotm_31 = 2.*quat_2*quat_4 - 2.*quat_1*quat_3;
    rotm_32 = 2.*quat_3*quat_4 + 2.*quat_1*quat_2;
    rotm_33 = quat_1*quat_1 - quat_2*quat_2 - quat_3*quat_3 + quat_4*quat_4;
  }

  int out_dim[] = {image_side, image_side, image_side};
  PyObject *out_array = (PyObject *)PyArray_FromDims(3, out_dim, NPY_FLOAT64);
  double *out = PyArray_DATA(out_array);

  //double edge_thickness = 1./image_side;
  double edge_thickness = 1.;
  double half_edge_thickness = edge_thickness/2.;
  const double phi = (1.+sqrt(5.))/2.; //golden ratio

  double corner1[] = {0., 1., phi};
  double corner2[] = {1., phi, 0.};
  double corner3[] = {phi, 0., 1.};

  double original_radius = sqrt(1. + pow(phi, 2));
  double size_scaling = radius/original_radius;

  corner1[0] *= size_scaling; corner1[1] *= size_scaling; corner1[2] *= size_scaling;
  corner2[0] *= size_scaling; corner2[1] *= size_scaling; corner2[2] *= size_scaling;
  corner3[0] *= size_scaling; corner3[1] *= size_scaling; corner3[2] *= size_scaling;
  
  double center_z = (corner1[0]+corner2[0]+corner3[0])/3.;
  double center_y = (corner1[1]+corner2[1]+corner3[1])/3.;
  double center_x = (corner1[2]+corner2[2]+corner3[2])/3.;

  double normal1_z = (corner1[0] + corner2[0])/2. - center_z;
  double normal1_y = (corner1[1] + corner2[1])/2. - center_y;
  double normal1_x = (corner1[2] + corner2[2])/2. - center_x;

  double normal2_z = (corner2[0] + corner3[0])/2. - center_z;
  double normal2_y = (corner2[1] + corner3[1])/2. - center_y;
  double normal2_x = (corner2[2] + corner3[2])/2. - center_x;

  double normal3_z = (corner3[0] + corner1[0])/2. - center_z;
  double normal3_y = (corner3[1] + corner1[1])/2. - center_y;
  double normal3_x = (corner3[2] + corner1[2])/2. - center_x;

  double edge_distance = sqrt(pow(normal1_z, 2) + pow(normal1_y, 2) + pow(normal1_x, 2));
  normal1_z /= edge_distance; normal1_y /= edge_distance; normal1_x /= edge_distance;
  normal2_z /= edge_distance; normal2_y /= edge_distance; normal2_x /= edge_distance;
  normal3_z /= edge_distance; normal3_y /= edge_distance; normal3_x /= edge_distance;

  double face_normal_3[] = {phi/3., 0., (2.*phi+1.)/3.};
  double face_normal_2[] = {(2.*phi+1.)/3., phi/3., 0.};
  double face_normal_1[] = {0., (2.*phi+1.)/3., phi/3.};
  double face_normal_center[] = {1., 1., 1.};

  double face_distance = sqrt(pow(center_z, 2) + pow(center_y, 2) + pow(center_x, 2))/size_scaling;
  face_normal_1[0] /= face_distance; face_normal_1[1] /= face_distance; face_normal_1[2] /= face_distance;
  face_normal_2[0] /= face_distance; face_normal_2[1] /= face_distance; face_normal_2[2] /= face_distance;
  face_normal_3[0] /= face_distance; face_normal_3[1] /= face_distance; face_normal_3[2] /= face_distance;
  face_normal_center[0] /= sqrt(3.); face_normal_center[1] /= sqrt(3.);
  face_normal_center[2] /= sqrt(3.);
  
  face_distance = sqrt(pow(center_z, 2) + pow(center_y, 2) + pow(center_x, 2));

  double no_rot_x, no_rot_y, no_rot_z;
  double x, y, z;
  double projected_x, projected_y, projected_z;
  double scalar_product, distance;
  double image_side_float = (double) image_side;
  int x_pixel, y_pixel, z_pixel, i;
  for (z_pixel = 0; z_pixel < image_side; z_pixel++) {
    //x = fabs(((double)x_pixel - image_side_float/2. + 0.5));
    no_rot_z = ((double)z_pixel - image_side_float/2. + 0.5);
    for (y_pixel = 0; y_pixel < image_side; y_pixel++) {
      //y = fabs(((double)y_pixel - image_side_float/2. + 0.5));
      no_rot_y = ((double)y_pixel - image_side_float/2. + 0.5);
      for (x_pixel = 0; x_pixel < image_side; x_pixel++) {
	//z = fabs(((double)z_pixel - image_side_float/2. + 0.5));
	no_rot_x = ((double)x_pixel - image_side_float/2. + 0.5);
	
	if (rotate == 1) {
	  //Transpose rotation matrix in comparison to Max' code since I don't rotate the icosahedron but the coordinate system.
	  z = fabs(no_rot_z*rotm_11 + no_rot_y*rotm_21 + no_rot_x*rotm_31);
	  y = fabs(no_rot_z*rotm_12 + no_rot_y*rotm_22 + no_rot_x*rotm_32);
	  x = fabs(no_rot_z*rotm_13 + no_rot_y*rotm_23 + no_rot_x*rotm_33);
	} else {
	  z = fabs(no_rot_z);
	  y = fabs(no_rot_y);
	  x = fabs(no_rot_x);
	}
	
	scalar_product = x*face_normal_center[2] + y*face_normal_center[1] + z*face_normal_center[0];
	projected_x = x * face_distance/scalar_product;
	projected_y = y * face_distance/scalar_product;
	projected_z = z * face_distance/scalar_product;

	i = z_pixel*image_side*image_side + y_pixel*image_side + x_pixel;
	if ((projected_x-center_x)*normal1_x + (projected_y-center_y)*normal1_y +
	    (projected_z-center_z)*normal1_z > edge_distance) {
	  distance = x*face_normal_1[2] + y*face_normal_1[1] + z*face_normal_1[0];

	  if (distance > face_distance + half_edge_thickness) {
	    out[i] = 0.;
	  } else if (distance < face_distance - half_edge_thickness) {
	    out[i] = 1.;
	  } else {
	    out[i] = 0.5 + (face_distance - distance) / edge_thickness;
	  }
	    
	} else if ((projected_x-center_x)*normal2_x + (projected_y-center_y)*normal2_y +
		   (projected_z-center_z)*normal2_z > edge_distance) {
	  distance = x*face_normal_2[2] + y*face_normal_2[1] + z*face_normal_2[0];
	  if (distance > face_distance + half_edge_thickness) {
	    out[i] = 0.;
	  } else if (distance < face_distance - half_edge_thickness) {
	    out[i] = 1.;
	  } else {
	    out[i] = 0.5 + (face_distance - distance) / edge_thickness;
	  }

	} else if ((projected_x-center_x)*normal3_x + (projected_y-center_y)*normal3_y +
		   (projected_z-center_z)*normal3_z > edge_distance) {
	  distance = x*face_normal_3[2] + y*face_normal_3[1] + z*face_normal_3[0];
	  if (distance > face_distance + half_edge_thickness) {
	    out[i] = 0.;
	  } else if (distance < face_distance - half_edge_thickness) {
	    out[i] = 1.;
	  } else {
	    out[i] = 0.5 + (face_distance - distance) / edge_thickness;
	  }

	} else {
	  distance = x*face_normal_center[2] + y*face_normal_center[1] + z*face_normal_center[0];
	  if (distance > face_distance + half_edge_thickness) {
	    out[i] = 0.;
	  } else if (distance < face_distance - half_edge_thickness) {
	    out[i] = 1.;
	  } else {
	    out[i] = 0.5 + (face_distance - distance) / edge_thickness;
	  }
	}
      }
    }
  }
  return out_array;
}

static PyMethodDef IcosahedronMethods[] = {
  {"icosahedron", (PyCFunction)icosahedron, METH_VARARGS|METH_KEYWORDS, icosahedron__doc__},
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

MOD_INIT(icosahedron)
{
  import_array();
  PyObject *m;
  MOD_DEF(m, "icosahedron", "Create icosahedron density map", IcosahedronMethods)
  if (m == NULL)
    return MOD_ERROR_VAL;
  return MOD_SUCCESS_VAL(m);
}
