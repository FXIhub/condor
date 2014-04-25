/*----------------------------------------------------------------------------------------------------- 
# PENGUIN 
# Simulator for diffractive single-particle imaging experiments with X-ray lasers
# http://xfel.icm.uu.se/penguin/
# ----------------------------------------------------------------------------------------------------- 
# Copyright 2014 Max Hantke, Filipe R.N.C. Maia, Tomas Ekeberg
# Penguin is distributed under the terms of the GNU General Public License
# ----------------------------------------------------------------------------------------------------- 
# General note:
#  All variables are in SI units by default. Exceptions explicit by variable name.
# ---------------------------------------------------------------------------------------------------*/

#include <Python.h>
#include <numpy/arrayobject.h>
#include <math.h>
#include <stdio.h>

PyDoc_STRVAR(icosahedron__doc__, "icosahedron(array_side, radius)\n\nGenerate an icosahedron. Radius is given in pixels and is defined as the distance to the corners.");
static PyObject *icosahedron(PyObject *self, PyObject *args, PyObject *kwargs)
{
  int image_side = 0.;
  double radius = 0.;

  if (!PyArg_ParseTuple(args, "id", &image_side, &radius)) {
    return NULL;
  }

  if (image_side <= 0) {
    PyErr_SetString(PyExc_ValueError, "Image side must be > 0.");
    return NULL;
  }

  if (radius <= 0.) {
    PyErr_SetString(PyExc_ValueError, "Radius must be > 0.");
  }

  int out_dim[] = {image_side, image_side, image_side};
  PyObject *out_array = (PyObject *)PyArray_FromDims(3, out_dim, NPY_FLOAT64);
  double *out = PyArray_DATA(out_array);

  //double edge_thickness = 4./(float)image_side;
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

  double x, y, z;
  double projected_x, projected_y, projected_z;
  double scalar_product, distance;
  double image_side_float = (double) image_side;
  for (int x_pixel = 0; x_pixel < image_side; x_pixel++) {
    x = fabs(((double)x_pixel - image_side_float/2. + 0.5));///image_side_float*4.);
    for (int y_pixel = 0; y_pixel < image_side; y_pixel++) {
      y = fabs(((double)y_pixel - image_side_float/2. + 0.5));///image_side_float*4.);
      for (int z_pixel = 0; z_pixel < image_side; z_pixel++) {
	z = fabs(((double)z_pixel - image_side_float/2. + 0.5));///image_side_float*4.);
	scalar_product = x*face_normal_center[2] + y*face_normal_center[1] + z*face_normal_center[0];
	projected_x = x * face_distance/scalar_product;
	projected_y = y * face_distance/scalar_product;
	projected_z = z * face_distance/scalar_product;

	if ((projected_x-center_x)*normal1_x + (projected_y-center_y)*normal1_y +
	    (projected_z-center_z)*normal1_z > edge_distance) {
	  distance = x*face_normal_1[2] + y*face_normal_1[1] + z*face_normal_1[0];
	  if (distance > face_distance + half_edge_thickness) {
	    out[x_pixel*image_side*image_side + y_pixel*image_side + z_pixel] = 0.;
	  } else if (distance < face_distance - half_edge_thickness) {
	    out[x_pixel*image_side*image_side + y_pixel*image_side + z_pixel] = 1.;
	  } else {
	    out[x_pixel*image_side*image_side + y_pixel*image_side + z_pixel] = 0.5 + (face_distance - distance) / edge_thickness;
	  }
	    
	} else if ((projected_x-center_x)*normal2_x + (projected_y-center_y)*normal2_y +
		   (projected_z-center_z)*normal2_z > edge_distance) {
	  distance = x*face_normal_2[2] + y*face_normal_2[1] + z*face_normal_2[0];
	  if (distance > face_distance + half_edge_thickness) {
	    out[x_pixel*image_side*image_side + y_pixel*image_side + z_pixel] = 0.;
	  } else if (distance < face_distance - half_edge_thickness) {
	    out[x_pixel*image_side*image_side + y_pixel*image_side + z_pixel] = 1.;
	  } else {
	    out[x_pixel*image_side*image_side + y_pixel*image_side + z_pixel] = 0.5 + (face_distance - distance) / edge_thickness;
	  }

	} else if ((projected_x-center_x)*normal3_x + (projected_y-center_y)*normal3_y +
		   (projected_z-center_z)*normal3_z > edge_distance) {
	  distance = x*face_normal_3[2] + y*face_normal_3[1] + z*face_normal_3[0];
	  if (distance > face_distance + half_edge_thickness) {
	    out[x_pixel*image_side*image_side + y_pixel*image_side + z_pixel] = 0.;
	  } else if (distance < face_distance - half_edge_thickness) {
	    out[x_pixel*image_side*image_side + y_pixel*image_side + z_pixel] = 1.;
	  } else {
	    out[x_pixel*image_side*image_side + y_pixel*image_side + z_pixel] = 0.5 + (face_distance - distance) / edge_thickness;
	  }

	} else {
	  distance = x*face_normal_center[2] + y*face_normal_center[1] + z*face_normal_center[0];
	  if (distance > face_distance + half_edge_thickness) {
	    out[x_pixel*image_side*image_side + y_pixel*image_side + z_pixel] = 0.;
	  } else if (distance < face_distance - half_edge_thickness) {
	    out[x_pixel*image_side*image_side + y_pixel*image_side + z_pixel] = 1.;
	  } else {
	    out[x_pixel*image_side*image_side + y_pixel*image_side + z_pixel] = 0.5 + (face_distance - distance) / edge_thickness;
	  }
	}
      }
    }
  }
  return out_array;
}

PyDoc_STRVAR(icosahedron_rot__doc__, "icosahedron(array_side, radius, rotation)\n\nGenerate an icosahedron. Radius is given in pixels and is defined as the distance to the corners. Rotation should be tuple of three euler angles");
static PyObject *icosahedron_rot(PyObject *self, PyObject *args, PyObject *kwargs)
{
  int image_side = 0.;
  double radius = 0.;
  PyObject *rotation_obj, *rotation_sequence;

  if (!PyArg_ParseTuple(args, "idO", &image_side, &radius, &rotation_obj)) {
    return NULL;
  }

  if (rotation_obj == NULL)
    return NULL;

  if (image_side <= 0) {
    PyErr_SetString(PyExc_ValueError, "Image side must be > 0.");
    return NULL;
  }

  if (radius <= 0.) {
    PyErr_SetString(PyExc_ValueError, "Radius must be > 0.");
  }

  rotation_sequence = PySequence_Fast(rotation_obj, "Expected a sequence");
  if (PySequence_Size(rotation_sequence) != 3) {
    PyErr_SetString(PyExc_ValueError, "Rotation must be of length 3");
    return NULL;
  }

  PyObject *seq_1 = PySequence_Fast_GET_ITEM(rotation_sequence, 0);
  PyObject *seq_2 = PySequence_Fast_GET_ITEM(rotation_sequence, 1);
  PyObject *seq_3 = PySequence_Fast_GET_ITEM(rotation_sequence, 2);
  const double euler_1 = PyFloat_AsDouble(seq_1);
  const double euler_2 = PyFloat_AsDouble(seq_2);
  const double euler_3 = PyFloat_AsDouble(seq_3);
  if (PyErr_Occurred()) {
    PyErr_SetString(PyExc_ValueError, "Rotation must be sequence of 3 euler angles (floats).");
    return NULL;
  }
    
  const double rotm_11 = cos(euler_2)*cos(euler_3);
  const double rotm_12 = -cos(euler_1)*sin(euler_3)+sin(euler_1)*sin(euler_2)*cos(euler_3);
  const double rotm_13 = sin(euler_1)*sin(euler_3)+cos(euler_1)*sin(euler_2)*cos(euler_3);
  const double rotm_21 = cos(euler_2)*sin(euler_3);
  const double rotm_22 = cos(euler_1)*cos(euler_3)+sin(euler_1)*sin(euler_2)*sin(euler_3);
  const double rotm_23 = -sin(euler_1)*cos(euler_3)+cos(euler_1)*sin(euler_2)*sin(euler_3);
  const double rotm_31 = -sin(euler_2);
  const double rotm_32 = sin(euler_1)*cos(euler_2);
  const double rotm_33 = cos(euler_1)*cos(euler_2);

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
  for (int x_pixel = 0; x_pixel < image_side; x_pixel++) {
    //x = fabs(((double)x_pixel - image_side_float/2. + 0.5));
    no_rot_x = ((double)x_pixel - image_side_float/2. + 0.5);
    for (int y_pixel = 0; y_pixel < image_side; y_pixel++) {
      //y = fabs(((double)y_pixel - image_side_float/2. + 0.5));
      no_rot_y = ((double)y_pixel - image_side_float/2. + 0.5);
      for (int z_pixel = 0; z_pixel < image_side; z_pixel++) {
	//z = fabs(((double)z_pixel - image_side_float/2. + 0.5));
	no_rot_z = ((double)z_pixel - image_side_float/2. + 0.5);

	//Transpose rotation matrix in comparison to Max' code since I don't rotate the icosahedron but the coordinate system.
	z = fabs(no_rot_z*rotm_11 + no_rot_y*rotm_21 + no_rot_x*rotm_31);
	y = fabs(no_rot_z*rotm_12 + no_rot_y*rotm_22 + no_rot_x*rotm_32);
	x = fabs(no_rot_z*rotm_13 + no_rot_y*rotm_23 + no_rot_x*rotm_33);
	
	scalar_product = x*face_normal_center[2] + y*face_normal_center[1] + z*face_normal_center[0];
	projected_x = x * face_distance/scalar_product;
	projected_y = y * face_distance/scalar_product;
	projected_z = z * face_distance/scalar_product;

	if ((projected_x-center_x)*normal1_x + (projected_y-center_y)*normal1_y +
	    (projected_z-center_z)*normal1_z > edge_distance) {
	  distance = x*face_normal_1[2] + y*face_normal_1[1] + z*face_normal_1[0];
	  if (distance > face_distance + half_edge_thickness) {
	    out[x_pixel*image_side*image_side + y_pixel*image_side + z_pixel] = 0.;
	  } else if (distance < face_distance - half_edge_thickness) {
	    out[x_pixel*image_side*image_side + y_pixel*image_side + z_pixel] = 1.;
	  } else {
	    out[x_pixel*image_side*image_side + y_pixel*image_side + z_pixel] = 0.5 + (face_distance - distance) / edge_thickness;
	  }
	    
	} else if ((projected_x-center_x)*normal2_x + (projected_y-center_y)*normal2_y +
		   (projected_z-center_z)*normal2_z > edge_distance) {
	  distance = x*face_normal_2[2] + y*face_normal_2[1] + z*face_normal_2[0];
	  if (distance > face_distance + half_edge_thickness) {
	    out[x_pixel*image_side*image_side + y_pixel*image_side + z_pixel] = 0.;
	  } else if (distance < face_distance - half_edge_thickness) {
	    out[x_pixel*image_side*image_side + y_pixel*image_side + z_pixel] = 1.;
	  } else {
	    out[x_pixel*image_side*image_side + y_pixel*image_side + z_pixel] = 0.5 + (face_distance - distance) / edge_thickness;
	  }

	} else if ((projected_x-center_x)*normal3_x + (projected_y-center_y)*normal3_y +
		   (projected_z-center_z)*normal3_z > edge_distance) {
	  distance = x*face_normal_3[2] + y*face_normal_3[1] + z*face_normal_3[0];
	  if (distance > face_distance + half_edge_thickness) {
	    out[x_pixel*image_side*image_side + y_pixel*image_side + z_pixel] = 0.;
	  } else if (distance < face_distance - half_edge_thickness) {
	    out[x_pixel*image_side*image_side + y_pixel*image_side + z_pixel] = 1.;
	  } else {
	    out[x_pixel*image_side*image_side + y_pixel*image_side + z_pixel] = 0.5 + (face_distance - distance) / edge_thickness;
	  }

	} else {
	  distance = x*face_normal_center[2] + y*face_normal_center[1] + z*face_normal_center[0];
	  if (distance > face_distance + half_edge_thickness) {
	    out[x_pixel*image_side*image_side + y_pixel*image_side + z_pixel] = 0.;
	  } else if (distance < face_distance - half_edge_thickness) {
	    out[x_pixel*image_side*image_side + y_pixel*image_side + z_pixel] = 1.;
	  } else {
	    out[x_pixel*image_side*image_side + y_pixel*image_side + z_pixel] = 0.5 + (face_distance - distance) / edge_thickness;
	  }
	}
      }
    }
  }
  return out_array;
}

static PyMethodDef IcosahedronMethods[] = {
  {"icosahedron", (PyCFunction)icosahedron_rot, METH_VARARGS, icosahedron_rot__doc__},
  //{"icosahedron_rot", (PyCFunction)icosahedron_rot, METH_VARARGS, icosahedron_rot__doc__},
  {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC initicosahedron(void)
{
  import_array();
  PyObject *m = Py_InitModule3("icosahedron", IcosahedronMethods, "Create icosahedron density map");
  if (m == NULL)
    return;
}
