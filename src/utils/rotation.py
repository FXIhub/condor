"""
This module is an implementation of a variety of tools for rotations in 3D space.
"""

import sys, numpy, types, pickle, time, math
 
import logging
logger = logging.getLogger(__name__)

from log import log_and_raise_error,log_warning,log_info,log_debug
import linalg

# CANONICAL ROTATION MATRICES   
# Rotation matrix around x-axis - observing the right hand rule
R_x = lambda t: numpy.array([[1., 0., 0.],
                             [0., numpy.cos(t), -numpy.sin(t)],
                             [0., numpy.sin(t), numpy.cos(t)]])
# Rotation matrix around y-axis - observing the right hand rule 
R_y = lambda t: numpy.array([[numpy.cos(t), 0., numpy.sin(t)],
                             [0., 1., 0.],
                             [-numpy.sin(t), 0., numpy.cos(t)]])
# Rotation matrix around z-axis - observing the right hand rule 
R_z = lambda t: numpy.array([[numpy.cos(t), -numpy.sin(t), 0.],
                             [numpy.sin(t), numpy.cos(t), 0.],
                             [0., 0., 1.]])

# Rotation of a given vector by a given angle with respect to one of the three principal axes
rot_x = lambda v,t: R_x(t).dot(v)
rot_y = lambda v,t: R_y(t).dot(v)
rot_z = lambda v,t: R_z(t).dot(v)

# CANONICAL QUATERNIONS
# Quaternion from angle and rotation unit vector coordinates (right-hand rule)
quat = lambda theta,ux,uy,uz: numpy.array([numpy.cos(theta/2.),
                                           numpy.sin(theta/2.)*ux,
                                           numpy.sin(theta/2.)*uy,
                                           numpy.sin(theta/2.)*uz])
# Quaternions for roations with respect to the x-, y- or z-axis
quat_x = lambda theta: quat(theta,1.,0.,0.)
quat_y = lambda theta: quat(theta,0.,1.,0.)
quat_z = lambda theta: quat(theta,0.,0.,1.)

class Rotation:
    r"""
    Class for a rotation in 3D space

    **Arguments:**

      :values (array): Array of values that define the rotation. For random rotations set values=``None`` and for example formalism=``'random'``. (default ``None``)

      :formalism: Formalism that defines how the argument values is interpreted. If ``None`` no rotation. (default ``None``)

        *Rotation formalism can be one of the following:*
    
        ======================== =========================================================================================================================== ===============================================================================
        ``formalism``            Variables                                                                                                                   ``values``
        ======================== =========================================================================================================================== ===============================================================================
        ``'quaternion'``         :math:`q = w + ix + jy + kz`                                                                                                :math:`[w,x,y,z]`
        ``'rotation_matrix'``    :math:`R = \begin{pmatrix} R_{11} & R_{12} & R_{13} \\ R_{21} & R_{22} & R_{23} \\ R_{31} & R_{32} & R_{33}  \end{pmatrix}` :math:`[[R_{11},R_{12},R_{13}],[R_{21},R_{22},R_{23}],[R_{31},R_{32},R_{33}]]`
        ``'euler_angles_zxz'``   :math:`e_1^{(z)}`, :math:`e_2^{(x)}`, :math:`e_3^{(z)}`                                                                     :math:`[e_1^{(y)},e_2^{(z)},e_3^{(y)}]`
        ``'euler_angles_xyx'``   :math:`e_1^{(x)}`, :math:`e_2^{(y)}`, :math:`e_3^{(x)}`                                                                     :math:`[e_1^{(x)},e_2^{(y)},e_3^{(x)}]`
        ``'euler_angles_xyz'``   :math:`e_1^{(x)}`, :math:`e_2^{(y)}`, :math:`e_3^{(z)}`                                                                     :math:`[e_1^{(x)},e_2^{(y)},e_3^{(z)}]`
        ``'euler_angles_yzx'``   :math:`e_1^{(y)}`, :math:`e_2^{(z)}`, :math:`e_3^{(x)}`                                                                     :math:`[e_1^{(y)},e_2^{(z)},e_3^{(x)}]`
        ``'euler_angles_zxy'``   :math:`e_1^{(z)}`, :math:`e_2^{(x)}`, :math:`e_3^{(y)}`                                                                     :math:`[e_1^{(z)},e_2^{(x)},e_3^{(y)}]`
        ``'euler_angles_zyx'``   :math:`e_1^{(z)}`, :math:`e_2^{(y)}`, :math:`e_3^{(x)}`                                                                     :math:`[e_1^{(z)},e_2^{(y)},e_3^{(x)}]`
        ``'euler_angles_yxz'``   :math:`e_1^{(y)}`, :math:`e_2^{(x)}`, :math:`e_3^{(z)}`                                                                     :math:`[e_1^{(y)},e_2^{(x)},e_3^{(z)}]`
        ``'euler_angles_xzy'``   :math:`e_1^{(x)}`, :math:`e_2^{(z)}`, :math:`e_3^{(y)}`                                                                     :math:`[e_1^{(x)},e_2^{(z)},e_3^{(y)}]`
        ``'random'``             *fully random rotation*                                                                                                     ``None``
        ``'random_x'``           *random rotation around* :math:`x` *axis*                                                                                   ``None``
        ``'random_y'``           *random rotation around* :math:`y` *axis*                                                                                   ``None``
        ``'random_z'``           *random rotation around* :math:`z` *axis*                                                                                   ``None``
        ======================== =========================================================================================================================== ===============================================================================
    """
    
    def __init__(self, values=None, formalism=None):
        self.rotation_matrix = None
        if values is None and formalism is None:
            # No rotation (rotation matrix = identity matrix)
            self.rotation_matrix = numpy.ones(shape=(3,3))
        elif formalism.startswith("euler_angles_") and len(formalism) == len("euler_angles_xyz"):
            self.set_with_euler_angles(values, rotation_axes=formalism[-3:])
        elif formalism == "rotation_matrix":
            self.set_with_rotation_matrix(values)
        elif formalism == "quaternion":
            self.set_with_quaternion(values)
        elif formalism in ["random","random_x","random_y","random_z"]:
            if values is not None:
                log_warning(logger, "Specified formalism=%s but values is not None." % formalism)

            self._set_as_random_formalism(formalism)
        else:
            log_and_raise_error(logger, "formalism=%s is not implemented" % formalism)
            return

    def set_with_euler_angles(self, euler_angles, rotation_axes="zxz"):
        r"""
        Set rotation with an array of three euler angles

        Args:
           :euler_angles: Array of the three euler angles representing consecutive rotations

        Kwargs:
           :rotation_axes(str): Rotation axes of the three consecutive rotations (default = \'zxz\') 
        """
        # Check input
        if euler_angles.size != 3:
            log_and_raise_error(logger, "Size of rotation variable does not match expected shape")
            return
        # Set rotation matrix
        self.rotation_matrix = euler_to_rotmx(euler_angles, rotation_axes)
        
    def set_with_rotation_matrix(self, rotation_matrix):
        r"""
        Set rotation with a rotation matrix

        Args:
           :rotation_matrix: 3x3 array representing the rotation matrix
        """
        # Check input
        if rotation_matrix.size != 9 or rotation_matrix.ndim != 2:
            log_and_raise_error(logger, "Size of rotation variable does not match expected shape")
            return
        I = rotation_matrix.dot(rotation_matrix.T)
        tol = 0.0001
        for i in range(3):
            if abs(I[i,i]-1.) > tol:
                log_and_raise_error(logger, "Given matrix cannot be a rotation matrix because it is not unitary")
        # Set rotation matrix
        self.rotation_matrix = rotation_matrix.copy()

    def set_with_quaternion(self, quaternion):
        r"""
        Set rotation with a quaternion

        Args:
           :quaternion: Numpy array representing the quaternion :math:`w+ix+jy+kz`: 

                        [:math:`w`, :math:`x`, :math:`y`, :math:`z`] = [:math:`\cos(\theta/2)`, :math:`u_x \sin(\theta/2)`, :math:`u_y \sin(\theta/2)`, :math:`u_z \sin(\theta/2)`] 

                        with :math:`\theta` being the rotation angle and :math:`\vec{u}=(u_x,u_y,u_z)` the unit vector that defines the axis of rotation.
        """
        # Check input
        if quaternion.size != 4:
            log_and_raise_error(logger, "Size of rotation variable does not match expected shape")
            return
        # Set rotation matrix
        self.rotation_matrix = rotmx_from_quat(quaternion)
        
    def _set_as_random_formalism(self, formalism):
        if formalism == "random":
            self.set_as_random()
        elif formalism == "random_x":
            self.set_as_random_x()
        elif formalism == "random_y":
            self.set_as_random_y()
        elif formalism == "random_z":
            self.set_as_random_z()
        
    def set_as_random(self):
        """
        Set new random rotation (fully random).
        """
        q = rand_quat()
        self.rotation_matrix = rotmx_from_quat(q)

    def set_as_random_x(self):
        """
        Set new random rotation around the :math:`x`-axis.
        """
        ang = numpy.random.rand()*2*numpy.pi
        self.rotation_matrix = R_x(ang)

    def set_as_random_y(self):
        """
        Set new random rotation around the :math:`y`-axis.
        """
        ang = numpy.random.rand()*2*numpy.pi
        self.rotation_matrix = R_y(ang)

    def set_as_random_z(self):
        """
        Set new random rotation around the :math:`z`-axis.
        """
        ang = numpy.random.rand()*2*numpy.pi
        self.rotation_matrix = R_z(ang)

    def invert(self):
        """
        Invert rotation
        """
        q = self.get_as_quaternion()
        q[1:] = -q[1:]
        self.set_with_quaternion(q)
        
    def is_similar(self, rotation, tol=0.00001):
        r"""
        Compare rotation with another instance of the Rotation class. If quaternion distance is smaller than tol return ``True``

        Args:
           :rotation (:class:`condor.utils.rotation.Rotation`): Instance of the Rotation class

        Kwargs:
           :tol (float): Tolerance for similarity. This is the maximum distance of the two quaternions in 4D space that will be interpreted for similar rotations. (default 0.00001)
        """
        q0 = self.get_as_quaternion(unique_representation=True)
        q1 = rotation.get_as_quaternion(unique_representation=True)
        err = numpy.sqrt(((q0-q1)**2).sum())
        return (err < tol)
            
    def rotate_vector(self, vector, order="xyz"):
        r"""
        Return the rotated copy of a given vector

        Args:
           :vector (array): 3D vector

        Kwargs:
           :order (str): Order of geometrical axes in array representation of the given vector (default ``'xyz'``)
        """
        # Check input
        if vector.size != 3 or vector.ndim != 1:
            log_and_raise_error(logger, "Cannot rotate vector. Vector has incompatible size (%i) or number of dimensions (%i)." % (vector.size,vector.ndim))
            return
        # Rotate
        if order == "xyz":
            return self.rotation_matrix.dot(vector)
        elif order == "zyx":
            return self.rotation_matrix.dot(vector[::-1])
        else:
            log_and_raise_error(logger, "Corrdinates in order=%s is invalid." % order)

    def rotate_vectors(self, vectors, order="xyz"):
        r"""
        Return the rotated copy of a given array of vectors

        Args:
           :vectors (array): Array of 3D vectors with shape (:math:`N`, 3) with :math:`N` denoting the number of 3D vectors

        Kwargs:
           :order (str): Order of geometrical axes in array representation of the given vector (default ``'xyz'``)
        """        
        # Check input
        if vectors.ndim != 2 and vectors.ndim != 1:
            log_and_raise_error(logger, "Cannot rotate vectors. Input must have either one or two dimensions")
            return           
        n_ax = list(vectors.shape)[-1]
        N = vectors.size
        Nv = N/3
        if vectors.ndim == 2 and n_ax != 3:
            log_and_raise_error(logger, "Cannot rotate vectors. The given array has length %i in last dimension but should be 3." % (n_ax))
            return
        if vectors.ndim == 1 and N % 3 != 0:
            log_and_raise_error(logger, "Cannot rotate vectors. The given array has size %i which is not a multiple of 3." % (n_ax))
            return
        # Rotate
        if order == "xyz":
            return numpy.array([numpy.dot(self.rotation_matrix,vectors.ravel()[i*3:(i+1)*3]) for i in numpy.arange(Nv)])
        elif order == "zyx":
            return numpy.array([numpy.dot(self.rotation_matrix,(vectors.ravel()[i*3:(i+1)*3])[::-1])[::-1] for i in numpy.arange(Nv)])
        else:
            log_and_raise_error(logger, "Corrdinates in order=%s is invalid." % order)

    def get_as_euler_angles(self, rotation_axes="zxz"):
        r"""
        Get rotation in Euler angle represantation :math:`[e_1^{(z)}, e_2^{(x)}, e_3^{(z)}]` (for the case of ``rotation_axis='zxz'``).

        Kwargs:
           :rotation_axes (str): Rotation axes of the three rotations (default ``'zxz'``) 
        """
        q = self.get_as_quaternion()
        return euler_from_quat(q, rotation_axes=rotation_axes)
    
    def get_as_rotation_matrix(self):
        r"""
        Get rotation in rotation matrix representation (3x3 array)
        """
        return self.rotation_matrix.copy()

    def get_as_quaternion(self, unique_representation=False):
        r"""
        Get rotation in quaternion representation :math:`[w, x, y, z]`.

        Kwargs:
           :unique_representation (bool): Make quaternion unique. For more details see the documentation of :func:`condor.utils.rotation.unique_representation_quat` (default = False)
        """
        q = quat_from_rotmx(self.rotation_matrix)
        if unique_representation:
            q = unique_representation_quat(q)
        return q

class Rotations:
    r"""
    Class for a list of rotations in 3D space

    Args:
      :values (array): Arrays of values that define the rotation. For random rotations set ``values = None`` (default ``None``)

      :formalism (str): See :class:`condor.utils.rotation.Rotation`. For no rotation set ``formalism = None`` (default ``None``)
    """    
    def __init__(self, values=None, formalism=None):
        """
       """
        if values is None and formalism is None:
            single = True
        elif formalism.startswith("euler_angles_") and len(formalism) == len("euler_angles_xyz"):
            values = numpy.asarray(values)
            single = values.ndim == 1
        elif formalism == "rotation_matrix":
            values = numpy.asarray(values)
            single = values.ndim == 2 
        elif formalism == "quaternion":
            values = numpy.asarray(values)
            single = values.ndim == 1
        elif formalism in ["random","random_x","random_y","random_z"]:
            single = True
        else:
            log_and_raise_error(logger, "formalism=%s is not implemented" % formalism)
            return
        self._formalism = formalism
        self._i = 0
        self._values = values
        # Initialise rotations
        if single:
            if values is None:
                # No rotation (rotation matrix = identity matrix)
                self._rotations = [Rotation()]
            else:
                self._rotations = [Rotation(values, formalism=formalism)]
        else:
            self._rotations = []
            for i in range(len(values)):
                self._rotations.append(Rotation(values[i], formalism=formalism))

    def get_formalism(self):
        """
        Return formalism that defines how the rotation values are geometrically interpreted
        """
        return self._formalism
                
    def get_next_rotation(self):
        """
        Iterate and return next rotation
        """
        if self._formalism in ["random","random_x","random_y","random_z"]:
            self._rotations[0]._set_as_random_formalism(self._formalism)
        rotation =  self.get_current_rotation()
        self._i += 1
        return rotation
    
    def get_current_rotation(self):
        """
        Return current rotation
        """
        return self._rotations[self._i % len(self._rotations)]

    def get_all_values(self):
        """
        Return all values that define the rotations
        """
        return self._values


# CONVERSIONS BETWEEN THE DIFFERENT REPRESENTATIONS

def euler_to_rotmx(euler_angles, rotation_axes="zxz"):
    r"""
    Obtain rotation matrix from three euler angles and the rotation axes

    Args:
      :euler_angles (array): Length-3 array of euler angles

    Kwargs:
      :rotation_axes (str): Rotation axes of the three consecutive Euler rotations (default ``'zxz'``) 
    """
    R = numpy.identity(3)
    for ang,ax in zip(euler_angles, rotation_axes):
        if ax == "x":
            R = R.dot(R_x(ang))
        elif ax == "y":
            R = R.dot(R_y(ang))
        elif ax == "z":
            R = R.dot(R_z(ang))
        else:
            log_and_raise_error(logger, "%s is not a valid axis" % ax)
            return
    return R

def rotmx_from_quat(q):
    r"""
    Create a rotation matrix from given quaternion ([Shoemake1992]_ page 128)

    Args:
       :quaternion (array): :math:`q = w + ix + jy + kz` (``values``: :math:`[w,x,y,z]`)

    The direction of rotation follows the right hand rule
    """
    w,x,y,z = q
    R = numpy.array([[1.-2.*(y**2+z**2),
                      2.*(x*y-w*z),
                      2.*(x*z+w*y)],
                     [2.*(x*y+w*z),
                      1.-2.*(x**2+z**2),
                      2.*(y*z-w*x)],
                     [2.*(x*z-w*y),
                      2.*(y*z+w*x),
                      1.-2.*(x**2+y**2)]])
    return R



# Quaternion from rotation matrix
def quat_from_rotmx(R):
    r"""
    Obtain the quaternion from a given rotation matrix (ref. [euclidianspace_mxToQuat]_)

    Args:
       :R: 3x3 array that represent the rotation matrix (see `Conventions <conventions.html#matrices>`_)
    """
    q = numpy.zeros(4, dtype="float")
    q[0] = numpy.sqrt( max( 0, 1 + R[0,0] + R[1,1] + R[2,2] ) ) / 2.
    q[1] = numpy.sqrt( max( 0, 1 + R[0,0] - R[1,1] - R[2,2] ) ) / 2.
    q[2] = numpy.sqrt( max( 0, 1 - R[0,0] + R[1,1] - R[2,2] ) ) / 2.
    q[3] = numpy.sqrt( max( 0, 1 - R[0,0] - R[1,1] + R[2,2] ) ) / 2.
    q[1] = numpy.copysign( q[1], R[2,1] - R[1,2] ) 
    q[2] = numpy.copysign( q[2], R[0,2] - R[2,0] ) 
    q[3] = numpy.copysign( q[3], R[1,0] - R[0,1] )
    return q

# Euler angles from rotation matrix
def euler_from_quat(q, rotation_axes="zxz"):
    r"""
    Return euler angles from quaternion (ref. [euclidianspace_mxToQuat]_, [euclidianspace_quatToEul]_).

    Args:
       :q: Numpy array :math:`[w,x,y,z]` that represents the quaternion

    Kwargs:
       :rotation_axes(str): Rotation axes of the three consecutive Euler rotations (default ``\'zxz\'``) 
    """
    if len(rotation_axes) != 3:
        print "Error: rotation_axes = %s is an invalid input." % rotation_axes
        return
    for s in rotation_axes:
        if s not in "xyz":
            print "Error: rotation_axes = %s is an invalid input." % rotation_axes
            return
    i1 = 0 if rotation_axes[0] == "x" else 1 if rotation_axes[0] == "y" else 2 if rotation_axes[0] == "z" else None
    i2 = 0 if rotation_axes[1] == "x" else 1 if rotation_axes[1] == "y" else 2 if rotation_axes[1] == "z" else None
    i3 = 0 if rotation_axes[2] == "x" else 1 if rotation_axes[2] == "y" else 2 if rotation_axes[2] == "z" else None
    v3 = numpy.array([0.,0.,0.])
    v3[i3] = 1.
    v3r = rotate_quat(v3, q)
    if ((i1==0) and (i2==2) and (i3==0)) or \
       ((i1==1) and (i2==0) and (i3==1)) or \
       ((i1==2) and (i2==1) and (i3==2)):
        e0 = numpy.arctan2(v3r[(i1+2)%3],v3r[(i1+1)%3])
        e1 = numpy.arccos(v3r[i1])
    elif ((i1==0) and (i2==2) and (i3==1)) or \
         ((i1==1) and (i2==0) and (i3==2)) or \
         ((i1==2) and (i2==1) and (i3==0)):
        e0 = numpy.arctan2(v3r[(i1+2)%3],v3r[(i1+1)%3])
        e1 = -numpy.arcsin(v3r[i1])
    elif ((i1==0) and (i2==1) and (i3==0)) or \
         ((i1==1) and (i2==2) and (i3==1)) or \
         ((i1==2) and (i2==0) and (i3==2)):
        e0 = numpy.arctan2(v3r[(i1+1)%3],-v3r[(i1+2)%3])
        e1 = numpy.arccos(v3r[i1])
    else:
        e0 = numpy.arctan2(-v3r[(i1+1)%3],v3r[(i1+2)%3])
        # The reference states this:
        #e1 = -numpy.arcsin(v3r[i1])
        # The tests only pass with the inverse sign, so I guess this is a typo.
        e1 = numpy.arcsin(v3r[i1])
    q1 = numpy.array([numpy.cos(e0/2.), 0., 0., 0.])
    q1[1+i1] = numpy.sin(e0/2.)
    q2 = numpy.array([numpy.cos(e1/2.), 0., 0., 0.])
    q2[1+i2] = numpy.sin(e1/2.)
    q12 = quat_mult(q1, q2)
    v3n = numpy.array([0.,0.,0.])
    v3n[(i3+1)%3] = 1.
    v3n12 = quat_vec_mult(q12, v3n)
    v3nG = quat_vec_mult(q, v3n)
    e2_mag = numpy.arccos(linalg.dotproduct(v3n12,v3nG))
    vc = linalg.crossproduct(v3n12, v3nG)
    m = linalg.dotproduct(vc, v3r)
    e2 = numpy.sign(m) * e2_mag
    return numpy.array([e0,e1,e2])



# CONVERSIONS FOR UNIQUE REPRESENTATION

def make_euler_unique_repax(euler):
    r"""
    Make Euler angles for roations with a repeated axis (zxz, zyz, yzy, yxy, xzx, xyx) unique (such that three euler angles uniquely define a rotation).

    *There are the following ambiguities:*

      1) Common circularity of rotations:

        .. math::

          (e_1,e_2,e_3) \Leftrightarrow (e_1 + n \cdot 2\pi ,e_2 + n \cdot 2\pi ,e_3 + n \cdot 2\pi) \, , \, n \in \mathbb{Z}

      2) Ambiguity resulting if first and last rotation is increased by :math:`\pi`:
   
        .. math::
    
          (e_1,e_2,e_3) \Leftrightarrow (e_1 + \pi, -e_2, e_3 + \pi)

      3) Gimbal lock if :math:`(e_1+e_3) \mod 2\pi = 0`:

        .. math::

          (e_1,e_2,e_3) \Leftrightarrow (a,e_2,-a)`, :math:`a \in \mathbb{Q}

    *The follwing conventions are being used (observe order!):*

      1) Allways:

        .. math::

          (e_1 \, \textrm{mod} \, 2\pi, e_2 \, \textrm{mod} \, 2\pi, e_3 \, \textrm{mod} \, 2\pi)

      2) If :math:`e_1 \geq \pi`:
 
        .. math::

          ((e_1 + pi) \mod 2\pi, 2\pi - e_2, (e_3 + \pi) \mod 2\pi)

      3) If :math:`| e_1+e_3-2\pi | < \epsilon`:
   
        .. math::

          (0,e_1,0)

    Args:
      euler (array): Numpy array with the three euler angles in radian.
    """
    pi = numpy.pi
    twopi = 2*pi
    euler = euler % twopi
    if euler[1] >= pi:
        euler[0] = (euler[0] + pi) % twopi
        euler[1] = twopi - euler[1]
        euler[2] = (euler[2] + pi) % twopi
    if abs(euler[0]+euler[2]-twopi) < numpy.finfo("float64").resolution:
        euler[0] = 0.
        euler[2] = 0.
    return euler

def unique_representation_quat(q):
    r"""
    Return the quaternion in a unique representation of rotations by avoiding the ambiguity that :math:`q` and :math:`-q` determine the same rotation.

    The convention is that the first non-zero coordinate value of the quaternion array has to be positive.

    Args:
      :q (array): Length-4 array [:math:`w`, :math:`x`, :math:`y`, :math:`z`] representing the qauternion :math:`q = w + ix + jy + kz`
    """
    if q[0] < 0:
        return -q
    elif q[0] == 0:
        if q[1] < 0:
            return -q
        elif q[1] == 0:
            if q[2] < 0:
                return -q
            elif q[2] == 0:
                if q[3] < 0:
                    return -q
    return q           

# QUATERNINON HELPER FUNCTIONS

# Normalisation of quaternion
def norm_quat(q, tolerance=0.00001):
    r"""
    Return a copy of the normalised quaternion (adjust length to 1 if deviation larger than given tolerance).

    Args:
       :tolerance(float): Maximum deviation of length before rescaling (default ``0.00001``)
    """
    # Adjust length
    l = linalg.length(q)
    if abs(l - 1.) > tolerance:
        q_norm = q_norm/linalg.length(q_norm)
    else:
        q_norm = q.copy()
    return q_norm


def quat_mult(q1, q2):
    r"""
    Return the product of two quaternions

    Args:
       :q1 (array): Length-4 array [:math:`w_1`, :math:`x_1`, :math:`y_1`, :math:`z_1`] that represents the first quaternion

       :q2 (array): Length-4 array [:math:`w_2`, :math:`x_2`, :math:`y_2`, :math:`z_2`] that represents the second quaternion
    """
    w1, x1, y1, z1 = q1
    w2, x2, y2, z2 = q2
    w = w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2
    x = w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2
    y = w1 * y2 + y1 * w2 + z1 * x2 - x1 * z2
    z = w1 * z2 + z1 * w2 + x1 * y2 - y1 * x2
    return numpy.array([w, x, y, z])

def quat_vec_mult(q, v):
    r"""
    Return the product of a quaternion and a vector

    Args:
       :q (array): Length-4 array [:math:`w`, :math:`x`, :math:`y`, :math:`z`] that represents the quaternion

       :v (array): Length-3 array [:math:`v_x`, :math:`v_y`, :math:`v_z`] that represents the vector
    """    
    q2 = numpy.array([0.,v[0],v[1],v[2]])
    return quat_mult(quat_mult(q, q2), quat_conj(q))[1:]

def quat_conj(q):
    r""" 
    Return the conjugate quaternion as a length-4 array [w,-ix,-jy,-kz]

    Args:
       :q (array): Numpy array :math:`[w,x,y,z]` that represents the quaternion
    """
    iq = q.copy()
    iq[1:] = -iq[1:]
    return iq

def rotate_quat(v, q):
    r"""
    Return rotated version of a given vector by a given quaternion

    Args:
       :v (array): Length-3 array :math:`[v_x,v_y,v_z]` that represents the vector      

       :q (array): Length-4 array [:math:`w`, :math:`x`, :math:`y`, :math:`z`] that represents the quaternion
    """
    return quat_vec_mult(q, v)

def rand_quat():
    r""" 
    Obtain a uniform random rotation in quaternion representation ([Shoemake1992]_ pages 129f)  
    """
    x0,x1,x2 = numpy.random.random(3)
    theta1 = 2.*numpy.pi*x1
    theta2 = 2.*numpy.pi*x2
    s1 = numpy.sin(theta1)
    s2 = numpy.sin(theta2)
    c1 = numpy.cos(theta1)
    c2 = numpy.cos(theta2)
    r1 = numpy.sqrt(1-x0)
    r2 = numpy.sqrt(x0)
    q = numpy.array([s1*r1, c1*r1, s2*r2, c2*r2])
    return q
