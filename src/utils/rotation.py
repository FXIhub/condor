import sys, numpy, types, pickle, time, math
 
import logging
logger = logging.getLogger(__name__)

from log import log_and_raise_error,log_warning,log_info,log_debug
import linalg

# Nomenclature:
# vector = [z,y,x]
# quaternion = [cos(theta/2), ux sin(theta/2), uy sin(theta/2), uz sin(theta/2)]
# We use right handed cartesian coordinates

    
class Rotation:
    
    def __init__(self, values=None, formalism=None):
        if values is None and formalism is None:
            # No rotation (rotation matrix = identity matrix)
            self.rotation_matrix = numpy.ones(shape=(3,3))
        elif formalism.startswith("euler_angles_") and len(formalism) == len("euler_angles_xyz"):
            self.set_with_euler_angles(values, axes=formalism[-3:])
        elif formalism == "rotation_matrix":
            self.set_with_rotation_matrix(values)
        elif formalism == "quaternion":
            self.set_with_quaternion(values)
        elif formalism in ["random","random_x","random_y","random_z"]:
            if values is not None:
                log_warning(logger, "Specified formalism=%s but values is not None." % formalism)
            self.set_random(formalism)
        else:
            log_and_raise_error(logger, "formalism=%s is not implemented" % formalism)
            return

    def set_with_euler_angles(self, euler_angles, axes="zxz"):
        # Check input
        if euler_angles.size != 3:
            log_and_raise_error(logger, "Size of rotation variable does not match expected shape")
            return
        # Set rotation matrix
        self.rotation_matrix = R_euler(euler_angles, axes)
        
    def set_with_rotation_matrix(self, rotation_matrix):
        # Check input
        if rotation_matrix.size != 9 or rotation_matrix.ndim != 2:
            log_and_raise_error(logger, "Size of rotation variable does not match expected shape")
            return
        I = rotation_matrix.dot(rotation_matrix.T)
        tol = 0.0001
        for Ii in I.ravel():
            if abs(Ii-1.) > tol:
                log_and_raise_error(logger, "Given matrix cannot be a rotation matrix because it is not unitary")
        # Set rotation matrix
        self.rotation_matrix = rotation_matrix.copy()

    def set_with_quaternion(self, quaternion):
        # Check input
        if quaternion.size != 4:
            log_and_raise_error(logger, "Size of rotation variable does not match expected shape")
            return
        # Set rotation matrix
        self.rotation_matrix = rotmx_from_quat(quaternion)

    def set_random(self):
        q = rand_quat()
        self.rotation_matrix = rotmx_from_quat(q)

    def set_random_x(self):
        ang = numpy.random.rand()*2*numpy.pi
        self.rotation_matrix = R_x(ang)

    def set_random_y(self):
        ang = numpy.random.rand()*2*numpy.pi
        self.rotation_matrix = R_y(ang)

    def set_random_z(self):
        ang = numpy.random.rand()*2*numpy.pi
        self.rotation_matrix = R_z(ang)

    def invert(self):
        q = self.get_as_quaternion()
        q[1:] = -q[1:]
        self.set_with_quaternion(q)

    def similar(self, rotation, tol=0.00001):
        q0 = self.get_as_quaternion(norm=True)
        q1 = rotation.get_as_quaternion(norm=True)
        err = numpy.sqrt(((q0-q1)**2).sum())
        return (err < tol)
            
    def rotate_vector(self, vector):
        # Check input
        if vector.size != 3 or vector.ndim != 1:
            log_and_raise_error(logger, "Cannot rotate vector. Vector has incompatible size (%i) or number of dimensions (%i)." % (vector.size,vector.ndim))
            return
        # Rotate
        return self.R.dot(vector)

    def rotate_vectors(self, vectors):
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
        return numpy.array([numpy.dot(self.rotation_matrix,vectors.ravel()[i*3:(i+1)*3]) for i in numpy.arange(Nv)])

    def get_as_euler_angles(self, axes="zxz"):
        q = self.get_as_quaternion()
        return euler_from_quat(q, axes=axes)
    
    def get_as_rotation_matrix(self):
        return self.rotation_matrix.copy()

    def get_as_quaternion(self, norm=True):
        q = quat_from_rotmx(self.rotation_matrix)
        if norm:
           q = norm_quat(q)
        return q
    

class Rotations:

    def __init__(self, values=None, formalism=None):
        if values is None and formalism is None:
            single = True
        elif formalism.startswith("euler_angles_") and len(formalism) == len("euler_angles_xyz"):
            single = values.ndim == 1
        elif formalism == "rotation_matrix":
            single = values.ndim == 2 
        elif formalism == "quaternion":
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
                self._rotations = [Rotation(values[i], formalism=formalism)]
        else:
            self._rotations = []
            for i in range(len(values)):
                self._rotations.append(Rotation(values[i], formalism=formalism))

    def get_formalism(self):
        return self._formalism
                
    def get_next(self):
        if self._formalism in ["random","random_x","random_y","random_z"]:
            if self._formalism == "random":
                self._rotations[0].set_random()
            elif self._formalism == "random_x":
                self._rotations[0].set_random_x()
            elif self._formalism == "random_y":
                self._rotations[0].set_random_y()
            elif self._formalism == "random_z":
                self._rotations[0].set_random_z()
        rotation =  self.get_current()
        self._i += 1
        return rotation
    
    def get_current(self):
        return self._rotations[self._i % len(self._rotations)]

    def get_all_values(self):
        return self._values
            
# Rotation matrix around x-axis observing the right hand rule
R_x = lambda t: numpy.array([[1., 0., 0.],
                             [0., numpy.cos(t), -numpy.sin(t)],
                             [0., numpy.sin(t), numpy.cos(t)]])
#R_x = lambda t: numpy.array([[numpy.cos(t), numpy.sin(t), 0.],
#                             [-numpy.sin(t), numpy.cos(t), 0.],
#                             [0., 0., 1.]])

# Rotation matrix around y-axis observing the right hand rule 
R_y = lambda t: numpy.array([[numpy.cos(t), 0., numpy.sin(t)],
                             [0., 1., 0.],
                             [-numpy.sin(t), 0., numpy.cos(t)]])
                             
#R_y = lambda t: numpy.array([[numpy.cos(t), 0., -numpy.sin(t)],
#                             [0., 1., 0.],
#                             [numpy.sin(t), 0., numpy.cos(t)]])

# Rotation matrix around z-axis observing the right hand rule 
R_z = lambda t: numpy.array([[numpy.cos(t), -numpy.sin(t), 0.],
                             [numpy.sin(t), numpy.cos(t), 0.],
                             [0., 0., 1.]])
#R_z = lambda t: numpy.array([[1., 0., 0.],
#                             [0., numpy.cos(t), numpy.sin(t)],
#                             [0., -numpy.sin(t), numpy.cos(t)]])

def R_euler(euler_angles, axes="zxz"):
    R = numpy.ones(shape=(3,3))
    for ang,ax in zip(euler_angles, axes):
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
            
# Extrinsic rotations:
# rotated vector = dot product of rotation matrix and vector
def rot_x(v,t):
    return R_x(t).dot(v)
def rot_y(v,t):
    return R_y(t).dot(v)
def rot_z(v,t):
    return R_z(t).dot(v)


# Uniform random rotations
# Reference: K. Shoemake. Uniform random rotations. In D. Kirk, editor, Graphics Gems III, pages 124-132. Academic, New York, 1992.
# Create a uniform random rotation quaternion (pages 129f)
def rand_quat():
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
# Create a rotation matrix from given quaternion (page 128)
def rotmx_from_quat(q):
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
    #R = numpy.array([[1.-2.*(x**2+y**2),
    #                  2.*(y*z+w*x),
    #                  2.*(x*z-w*y)],
    #                 [2.*(y*z-w*x),
    #                  1.-2.*(x**2+z**2),
    #                  2.*(x*y+w*z)],
    #                 [2.*(x*z+w*y),
    #                  2.*(x*y-w*z),
    #                  1.-2.*(y**2+z**2)]])
    return R

# Quaternion from angle and rotation unit vector coordinates
quat = lambda theta,ux,uy,uz: numpy.array([numpy.cos(theta/2.),
                                           numpy.sin(theta/2.)*ux,
                                           numpy.sin(theta/2.)*uy,
                                           numpy.sin(theta/2.)*uz])
# Quaternions for roations with respect to the x-, y- or z-axis
quatx = lambda theta: quat(theta,1.,0.,0.)
quaty = lambda theta: quat(theta,0.,1.,0.)
quatz = lambda theta: quat(theta,0.,0.,1.)

# Normalisation of quaternion that avoids ambiguous quaternion roations (in agreement with Tomas' EMC code)
def norm_quat(q, tolerance=0.00001):
    # Length
    q_norm = q.copy()
    l = linalg.length(q)
    if abs(l - 1.) > tolerance:
        q_norm = q_norm/linalg.length(q_norm)
    else:
        q_norm = q
    # Signs
    if q_norm[0] < 0:
        q_norm = -q_norm
    elif q_norm[0] == 0:
        if q_norm[1] < 0:
            q_norm = -q_norm
        elif q_norm[1] == 0:
            if q_norm[2] < 0:
                q_norm = -q_norm
            elif q_norm[2] == 0:
                if q_norm[3] < 0:
                    q_norm = -q_norm
    return q_norm           
# Test
def _test_rotmx_from_quat(N=1000):
    err = numpy.ones(N)
    for i in range(N):
        theta = numpy.random.rand()*2*numpy.pi
        #print "R_x"
        R_x0 = R_x(theta)
        R_x1 = rotmx_from_quat(quatx(theta))
        errx = ((R_x0-R_x1)**2).sum()
        #print "Error = %e" % errx
        #print "R_y"
        R_y0 = R_y(theta)
        R_y1 = rotmx_from_quat(quaty(theta))
        erry = ((R_y0-R_y1)**2).sum()
        #print "Error = %e" % erry
        #print "R_z"
        R_z0 = R_z(theta)
        R_z1 = rotmx_from_quat(quatz(theta))
        errz = ((R_z0-R_z1)**2).sum()
        err[i] = errx+erry+errz
    if err.max() < 0.00001:
        print "Test successful (max. error = %e with %i samples)." % (err.max(),N)
    else:
        print "Test failed (max. error = %e with %i samples)." % (err.max(),N)

def inv_quat(quat):
    iquat = quat.copy()
    iquat[1:] = -iquat[1:]
    return iquat

def quat_mult(q1, q2):
    w1, x1, y1, z1 = q1
    w2, x2, y2, z2 = q2
    w = w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2
    x = w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2
    y = w1 * y2 + y1 * w2 + z1 * x2 - x1 * z2
    z = w1 * z2 + z1 * w2 + x1 * y2 - y1 * x2
    return numpy.array([w, x, y, z])

def quat_vec_mult(q1, v1):
    q2 = numpy.array([0.,v1[2],v1[1],v1[0]])
    return quat_mult(quat_mult(q1, q2), quat_conj(q1))[1:][::-1]

def quat_conj(q):
    w, x, y, z = q
    return (w, -x, -y, -z)

def rotate_quat(v,q):
    return quat_vec_mult(q, v)
    
# Quaternion from rotation matrix
# Method from http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/
def quat_from_rotmx(R):
    q = numpy.zeros(4, dtype="float")
    q[0] = numpy.sqrt( max( 0, 1 + R[0,0] + R[1,1] + R[2,2] ) ) / 2.
    q[1] = numpy.sqrt( max( 0, 1 + R[0,0] - R[1,1] - R[2,2] ) ) / 2.
    q[2] = numpy.sqrt( max( 0, 1 - R[0,0] + R[1,1] - R[2,2] ) ) / 2.
    q[3] = numpy.sqrt( max( 0, 1 - R[0,0] - R[1,1] + R[2,2] ) ) / 2.
    q[1] = numpy.copysign( q[1], R[2,1] - R[1,2] ) 
    q[2] = numpy.copysign( q[2], R[0,2] - R[2,0] ) 
    q[3] = numpy.copysign( q[3], R[1,0] - R[0,1] )
    return q
# Test
def _test_quat_from_rotmx(N=1000):
    err = numpy.ones(N)
    for i in range(N):
        q0 = rand_quat()
        q0 = norm_quat(q0)
        R = rotmx_from_quat(q0)
        q1 = quat_from_rotmx(R)
        err[i] = abs(q0-q1).sum()
    if err.max() < 0.00001:
        print "Test successful (max. error = %e with %i samples)." % (err.max(),N)
    else:
        print "Test failed (max. error = %e with %i samples)." % (err.max(),N)

# Euler angles from rotation matrix
# Sources:
# http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/
# http://www.euclideanspace.com/maths/geometry/rotations/conversions/quaternionToEuler/quat_2_euler_paper_ver2-1.pdf  
def euler_from_quat(q, axes="zxz"):
    if len(axes) != 3:
        print "Error: axes = %s is an invalid input." % axes
        return
    for s in axes:
        if s not in "xyz":
            print "Error: axes = %s is an invalid input." % axes
            return
    i1 = 2 if axes[0] == "x" else 1 if axes[0] == "y" else 0 if axes[0] == "z" else None
    i2 = 2 if axes[1] == "x" else 1 if axes[1] == "y" else 0 if axes[1] == "z" else None
    i3 = 2 if axes[2] == "x" else 1 if axes[2] == "y" else 0 if axes[2] == "z" else None
    v3 = numpy.array([0.,0.,0.])
    v3[i3] = 1.
    v3r = rotate_quat(v3, q)
    if ((i1==2) and (i2==0) and (i3==2)) or \
       ((i1==1) and (i2==2) and (i3==1)) or \
       ((i1==0) and (i2==1) and (i3==0)):
        e0 = numpy.arctan2(v3r[(i1-2)%3],v3r[(i1-1)%3])
        e1 = numpy.arccos(v3r[i1])
    elif ((i1==2) and (i2==0) and (i3==1)) or \
         ((i1==1) and (i2==2) and (i3==0)) or \
         ((i1==0) and (i2==1) and (i3==2)):
        e0 = numpy.arctan2(v3r[(i1-2)%3],v3r[(i1-1)%3])
        e1 = -numpy.arcsin(v3r[i1])
    elif ((i1==2) and (i2==1) and (i3==2)) or \
         ((i1==1) and (i2==0) and (i3==1)) or \
         ((i1==0) and (i2==2) and (i3==0)):
        e0 = numpy.arctan2(v3r[(i1-1)%3],-v3r[(i1-2)%3])
        e1 = numpy.arccos(v3r[i1])
    else:
        #e0 = numpy.arctan2(-v3r[(i1-1)%3],v3r[(i1-2)%3])
        #e1 = -numpy.arcsin(v3r[i1])
        e0 = numpy.arctan2(-v3r[(i1-1)%3],v3r[(i1-2)%3])
        e1 = numpy.arcsin(v3r[i1])
    e02 = numpy.arctan2(v3r[2],-v3r[1])
    e12 = numpy.arccos(v3r[0])
    q1 = numpy.array([numpy.cos(e0/2.), 0., 0., 0.])
    q1[3-i1] = numpy.sin(e0/2.)
    q2 = numpy.array([numpy.cos(e1/2.), 0., 0., 0.])
    q2[3-i2] = numpy.sin(e1/2.)
    q12 = quat_mult(q1, q2)
    v3n = numpy.array([0.,0.,0.])
    v3n[(i3-1)%3] = 1.
    v3n12 = quat_vec_mult(q12, v3n)
    v3nG = quat_vec_mult(q, v3n)
    e2_mag = numpy.arccos(linalg.dotproduct(v3n12,v3nG))
    vc = linalg.crossproduct(v3n12[::-1], v3nG[::-1])[::-1]
    m = linalg.dotproduct(vc, v3r)
    e2 = numpy.sign(m) * e2_mag
    return numpy.array([e0,e1,e2])

def norm_euler_repax(euler):
    """
    Normalising Euler angles for roations with a repeated axis (zxz, zyz, yzy, yxy, xzx, xyx) such that a set uniquely defines a rotation.
    There are the following ambiguities:
    1) (euler[0],euler[1],euler[2]) <=> (euler[0]+n*2*pi,euler[1]+n*2*pi,euler[2]+n*2*pi) , n: any integer
    2) (euler[0],euler[1],euler[2]) <=> (euler[0]+pi,-euler[1],euler[2]+pi)
    3) If ((euler[0]+euler[2]) mod 2pi) == 0: (euler[0],euler[1],euler[2]) <=> (a,euler[1],-a), a rational 

    The follwing conventions are being used (observe order!):
    1) Modulo 2pi of all angles
    2) If euler[0] >= pi: 
       a) euler[0] = (euler[0] + pi) modulo 2pi
       b) euler[1] = 2pi - euler[1]
       c) euler[2] = (euker[2] + pi) modulo 2pi
    3) If | euler[0]+euler[2]-2pi | < epsilon:
       euler[0] = 0.
       euler[2] = 0.
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
        
def _test_euler_from_quat(N=1000):
    axes = ["xyz",
            "xzy",
            "yxz",
            "yzx",
            "zxy",
            "zyx",
            "yxy",
            "zxz",
            "xyx",
            "zyz",
            "xzx",
            "yzy"]
    succ = 0
    NN = len(axes)*N
    err = numpy.ones(NN)
    k = 0
    for ax in axes:
        for i in range(N):
            v0 = numpy.random.random(3)
            euler0 = numpy.random.random(3) * numpy.pi * 2
            euler0[1] /= 2.
            #euler0 = norm_euler_zxz(euler0)
            exec "R0 = R_%s(euler0[0]).dot(R_%s(euler0[1]).dot(R_%s(euler0[2])))" % (ax[0],ax[1],ax[2])
            # For testing
            q0 = quat_from_rotmx(R0)
            q0 = norm_quat(q0)
            # Important call
            euler1 = euler_from_quat(q0, ax)
            # ----
            #euler1 = norm_euler_zxz(euler1)
            exec "R1 = R_%s(euler1[0]).dot(R_%s(euler1[1]).dot(R_%s(euler1[2])))" % (ax[0],ax[1],ax[2])
            q1 = quat_from_rotmx(R1)
            q1 = norm_quat(q1)
            #de = euler0-euler1
            dq = q0-q1
            #err[k] = abs(de).sum()
            err[k] = abs(dq).sum()
            if err[k] > 0.00001:
                print "ERROR in ax %s" % ax
            else:
                succ += 1
            k += 1
    if succ == NN:
        print "Test successful (max. error = %e, succeeded in %i/%i cases)." % (err.max(), succ, NN)
    else:
        print "Test failed (max. error = %e, succeeded in %i/%i cases)." % (err.max(), succ, NN)

# Rotation matrix Rxyz = Rx Ry Rz
# If a vector is pre-multiplied by Rxyz the operation represents
# - consecutive intrinsic roations x-y'-z'' OR
# - consecutive extrinsic rotations z-y-x.
# The rotation angles of the three consecutive rotations are E0, E1, E2 repectively. 
#R_xyz = lambda E0,E1,E2: numpy.array([[numpy.cos(E0)*numpy.cos(E1),
#                                      numpy.cos(E2)*numpy.sin(E0)+numpy.cos(E0)*numpy.sin(E1)*numpy.sin(E2),
#                                      numpy.sin(E0)*numpy.sin(E2)-numpy.cos(E0)*numpy.cos(E2)*numpy.sin(E1)],
#                                     [-numpy.cos(E1)*numpy.sin(E0),
#                                      numpy.cos(E0)*numpy.cos(E2)-numpy.sin(E0)*numpy.sin(E1)*numpy.sin(E2),
#                                      numpy.cos(E0)*numpy.sin(E2)+numpy.cos(E2)*numpy.sin(E0)*numpy.sin(E1)],
#                                     [numpy.sin(E1),
#                                      -numpy.cos(E1)*numpy.sin(E2),
#                                      numpy.cos(E1)*numpy.cos(E2)]])


# Rotation matrix Rzxz = Rz Rx Rz
# If a vector is pre-multiplied by Rzxz the operation represents
# - consecutive intrinsic roations z-x'-z'' OR
# - consecutive extrinsic rotations z-x-z.
# The rotation angles of the three consecutive rotations are E0, E1, E2 repectively. 
#Rzxz = lambda E0,E1,E2: numpy.array([[numpy.sin(E0)*numpy.sin(E1),
#                                      -numpy.cos(E0)*numpy.sin(E2)-numpy.sin(E0)*numpy.cos(E1)*numpy.cos(E2),
#                                      numpy.cos(E0)*numpy.cos(E2)-numpy.sin(E0)*numpy.cos(E1)*numpy.sin(E2)],
#                                     [-numpy.cos(E0)*numpy.sin(E1),
#                                      -numpy.sin(E0)*numpy.sin(E2)+numpy.cos(E0)*numpy.cos(E1)*numpy.cos(E2),
#                                      numpy.sin(E0)*numpy.cos(E0)+numpy.cos(E0)*numpy.cos(E1)*numpy.sin(E2)],
#                                     [numpy.sin(E0)*numpy.sin(E1),
#                                      -numpy.cos(E0)*numpy.sin(E2)-numpy.sin(E0)*numpy.cos(E1)*numpy.cos(E2),
#                                      numpy.cos(E0)*numpy.cos(E2)-numpy.sin(E0)*numpy.cos(E1)*numpy.sin(E2)]])
#R_zxz = lambda E0,E1,E2: numpy.array([[numpy.cos(E1),
#                                      numpy.sin(E1)*numpy.cos(E2),
#                                      numpy.sin(E1)*numpy.sin(E2)],
#                                     [-numpy.cos(E0)*numpy.sin(E1),
#                                      -numpy.sin(E0)*numpy.sin(E2)+numpy.cos(E0)*numpy.cos(E1)*numpy.cos(E2),
#                                      numpy.sin(E0)*numpy.cos(E2)+numpy.cos(E0)*numpy.cos(E1)*numpy.sin(E2)],
#                                     [numpy.sin(E0)*numpy.sin(E1),
#                                      -numpy.cos(E0)*numpy.sin(E2)-numpy.sin(E0)*numpy.cos(E1)*numpy.cos(E2),
#                                      numpy.cos(E0)*numpy.cos(E2)-numpy.sin(E0)*numpy.cos(E1)*numpy.sin(E2)]])
                                     


    
#def rotation(v,e0,e1,e2,convention="zxz"):
#    if convention == "zxz":
#        R = R_zxz
#    elif convention == "xyz":
#        R = R_xyz
#    else:
#        print "ERROR: Rotation convention=%s is not implemented." % convention
#        return
#    return R(e0,e1,e2).dot(v)      
    

# Rotation in libspimage
#def R_old(E0,E1,E2):
#    R_old = numpy.ones(shape=(3,3))
#    R_old[ 0, 0] = numpy.cos(E1)*numpy.cos(E2)
#    R_old[ 0, 1] = -numpy.cos(E0)*numpy.sin(E2)+numpy.sin(E0)*numpy.sin(E1)*numpy.cos(E2)
#    R_old[ 0, 2] =  numpy.sin(E0)*numpy.sin(E2)+numpy.cos(E0)*numpy.sin(E1)*numpy.cos(E2)
#    R_old[ 1, 0] = numpy.cos(E1)*numpy.sin(E2)
#    R_old[ 1, 1] = numpy.cos(E0)*numpy.cos(E2)+numpy.sin(E0)*numpy.sin(E1)*numpy.sin(E2)
#    R_old[ 1, 2] = -numpy.sin(E0)*numpy.cos(E2)+numpy.cos(E0)*numpy.sin(E1)*numpy.sin(E2)
#    R_old[ 2, 0] = -numpy.sin(E1)
#    R_old[ 2, 1] = numpy.sin(E0)*numpy.cos(E1)
#    R_old[ 2, 2] = numpy.cos(E0)*numpy.cos(E1)
#    return R_old

#def random_euler_angles():
#    """
#    Generates a triplet (phi, theta, psi) of random Euler angles.
#    """
#    r1,r2,r3 = numpy.random.random(3)
#    q1 = numpy.sqrt(1.0-r1)*numpy.sin(2.0*numpy.pi*r2)
#    q2 = numpy.sqrt(1.0-r1)*numpy.cos(2.0*numpy.pi*r2)
#    q3 = numpy.sqrt(r1)*numpy.sin(2.0*numpy.pi*r3)
#    q4 = numpy.sqrt(r1)*numpy.cos(2.0*numpy.pi*r3)
#    e1 = math.atan2(2.0*(q1*q2+q3*q4), 1.0-2.0*(q2**2+q3**2))
#    e2 = math.asin(2.0*(q1*q3-q4*q2))
#    e3 = math.atan2(2.0*(q1*q4+q2*q3), 1.0-2.0*(q3**2+q4**2))
#x    return (e1,e2,e3)

