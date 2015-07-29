# -----------------------------------------------------------------------------------------------------
# CONDOR
# Simulator for diffractive single-particle imaging experiments with X-ray lasers
# http://xfel.icm.uu.se/condor/
# -----------------------------------------------------------------------------------------------------
# Copyright 2014 Max Hantke, Filipe R.N.C. Maia, Tomas Ekeberg
# Condor is distributed under the terms of the GNU General Public License
# -----------------------------------------------------------------------------------------------------
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but without any warranty; without even the implied warranty of
# merchantability or fitness for a pariticular purpose. See the
# GNU General Public License for more details.
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
# -----------------------------------------------------------------------------------------------------
# General note:
# All variables are in SI units by default. Exceptions explicit by variable name.
# -----------------------------------------------------------------------------------------------------

import sys, numpy, types, pickle, time, math
 
import logging
logger = logging.getLogger("Condor")
from log import log

# Nomenclature:
# vector = [z,y,x]
# We use right handed cartesian coordinates


# Rotation matrix around x-axis observing the right hand rule
Rx = lambda t: numpy.array([[numpy.cos(t), numpy.sin(t), 0.],
                            [-numpy.sin(t), numpy.cos(t), 0.],
                            [0., 0., 1.]])

# Rotation matrix around y-axis observing the right hand rule 
Ry = lambda t: numpy.array([[numpy.cos(t), 0., -numpy.sin(t)],
                            [0., 1., 0.],
                            [numpy.sin(t), 0., numpy.cos(t)]])

# Rotation matrix around z-axis observing the right hand rule 
Rz = lambda t: numpy.array([[1., 0., 0.],
                            [0., numpy.cos(t), numpy.sin(t)],
                            [0., -numpy.sin(t), numpy.cos(t)]])

# Extrinsic rotations:
# rotated vector = dot product of rotation matrix and vector
def rot_x(v,t):
    return Rx(t).dot(v)
def rot_y(v,t):
    return Ry(t).dot(v)
def rot_z(v,t):
    return Rz(t).dot(v)


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
    R = numpy.array([[1.-2.*(x**2+y**2),
                      2.*(y*z+w*x),
                      2.*(x*z-w*y)],
                     [2.*(y*z-w*x),
                      1.-2.*(x**2+z**2),
                      2.*(x*y+w*z)],
                     [2.*(x*z+w*y),
                      2.*(x*y-w*z),
                      1.-2.*(y**2+z**2)]])
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

def normalize(v, tolerance=0.00001):
    mag2 = sum(n * n for n in v)
    if abs(mag2 - 1.0) > tolerance:
        mag = sqrt(mag2)
        v = tuple(n / mag for n in v)
    return v
def norm_quat(q, tolerance=0.00001):
    # Length
    q_norm = q.copy()
    l = length(q)
    if abs(l - 1.) > tolerance:
        q_norm = q_norm/length(q_norm)
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
        #print "Rx"
        Rx0 = Rx(theta)
        Rx1 = rotmx_from_quat(quatx(theta))
        errx = ((Rx0-Rx1)**2).sum()
        #print "Error = %e" % errx
        #print "Ry"
        Ry0 = Ry(theta)
        Ry1 = rotmx_from_quat(quaty(theta))
        erry = ((Ry0-Ry1)**2).sum()
        #print "Error = %e" % erry
        #print "Rz"
        Rz0 = Rz(theta)
        Rz1 = rotmx_from_quat(quatz(theta))
        errz = ((Rz0-Rz1)**2).sum()
        err[i] = errx+erry+errz
    if err.max() < 0.00001:
        print "Test successful (max. error = %e with %i samples)." % (err.max(),N)
    else:
        print "Test failed (max. error = %e with %i samples)." % (err.max(),N)

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
    q[0] = numpy.sqrt( max( 0, 1 + R[2,2] + R[1,1] + R[0,0] ) ) / 2.
    q[1] = numpy.sqrt( max( 0, 1 + R[2,2] - R[1,1] - R[0,0] ) ) / 2.
    q[2] = numpy.sqrt( max( 0, 1 - R[2,2] + R[1,1] - R[0,0] ) ) / 2.
    q[3] = numpy.sqrt( max( 0, 1 - R[2,2] - R[1,1] + R[0,0] ) ) / 2.
    q[1] = numpy.copysign( q[1], R[0,1] - R[1,0] ) 
    q[2] = numpy.copysign( q[2], R[2,0] - R[0,2] ) 
    q[3] = numpy.copysign( q[3], R[1,2] - R[2,1] )
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
def euler_from_quat(q, mode="zxz"):
    if len(mode) != 3:
        print "Error: mode = %s is an invalid input." % mode
        return
    for s in mode:
        if s not in "xyz":
            print "Error: mode = %s is an invalid input." % mode
            return
    i1 = 2 if mode[0] == "x" else 1 if mode[0] == "y" else 0 if mode[0] == "z" else None
    i2 = 2 if mode[1] == "x" else 1 if mode[1] == "y" else 0 if mode[1] == "z" else None
    i3 = 2 if mode[2] == "x" else 1 if mode[2] == "y" else 0 if mode[2] == "z" else None
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
    e2_mag = numpy.arccos(dotproduct(v3n12,v3nG))
    vc = crossproduct(v3n12[::-1], v3nG[::-1])[::-1]
    m = dotproduct(vc, v3r)
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
    modes = ["xyz",
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
    NN = len(modes)*N
    err = numpy.ones(NN)
    k = 0
    for mode in modes:
        for i in range(N):
            v0 = numpy.random.random(3)
            euler0 = numpy.random.random(3) * numpy.pi * 2
            euler0[1] /= 2.
            #euler0 = norm_euler_zxz(euler0)
            exec "R0 = R%s(euler0[0]).dot(R%s(euler0[1]).dot(R%s(euler0[2])))" % (mode[0],mode[1],mode[2])
            # For testing
            q0 = quat_from_rotmx(R0)
            q0 = norm_quat(q0)
            # Important call
            euler1 = euler_from_quat(q0, mode)
            # ----
            #euler1 = norm_euler_zxz(euler1)
            exec "R1 = R%s(euler1[0]).dot(R%s(euler1[1]).dot(R%s(euler1[2])))" % (mode[0],mode[1],mode[2])
            q1 = quat_from_rotmx(R1)
            q1 = norm_quat(q1)
            #de = euler0-euler1
            dq = q0-q1
            #err[k] = abs(de).sum()
            err[k] = abs(dq).sum()
            if err[k] > 0.00001:
                print "ERROR in mode %s" % mode
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
Rxyz = lambda E0,E1,E2: numpy.array([[numpy.cos(E0)*numpy.cos(E1),
                                      numpy.cos(E2)*numpy.sin(E0)+numpy.cos(E0)*numpy.sin(E1)*numpy.sin(E2),
                                      numpy.sin(E0)*numpy.sin(E2)-numpy.cos(E0)*numpy.cos(E2)*numpy.sin(E1)],
                                     [-numpy.cos(E1)*numpy.sin(E0),
                                      numpy.cos(E0)*numpy.cos(E2)-numpy.sin(E0)*numpy.sin(E1)*numpy.sin(E2),
                                      numpy.cos(E0)*numpy.sin(E2)+numpy.cos(E2)*numpy.sin(E0)*numpy.sin(E1)],
                                     [numpy.sin(E1),
                                      -numpy.cos(E1)*numpy.sin(E2),
                                      numpy.cos(E1)*numpy.cos(E2)]])


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
Rzxz = lambda E0,E1,E2: numpy.array([[numpy.cos(E1),
                                      numpy.sin(E1)*numpy.cos(E2),
                                      numpy.sin(E1)*numpy.sin(E2)],
                                     [-numpy.cos(E0)*numpy.sin(E1),
                                      -numpy.sin(E0)*numpy.sin(E2)+numpy.cos(E0)*numpy.cos(E1)*numpy.cos(E2),
                                      numpy.sin(E0)*numpy.cos(E2)+numpy.cos(E0)*numpy.cos(E1)*numpy.sin(E2)],
                                     [numpy.sin(E0)*numpy.sin(E1),
                                      -numpy.cos(E0)*numpy.sin(E2)-numpy.sin(E0)*numpy.cos(E1)*numpy.cos(E2),
                                      numpy.cos(E0)*numpy.cos(E2)-numpy.sin(E0)*numpy.cos(E1)*numpy.sin(E2)]])
                                     


    
def rotation(v,e0,e1,e2,convention="zxz"):
    if convention == "zxz":
        R = Rzxz
    elif convention == "xyz":
        R = Rxyz
    else:
        print "ERROR: Rotation convention=%s is not implemented." % convention
        return
    return R(e0,e1,e2).dot(v)      
    

# Rotation in libspimage
def R_old(E0,E1,E2):
    R_old = numpy.ones(shape=(3,3))
    R_old[ 0, 0] = numpy.cos(E1)*numpy.cos(E2)
    R_old[ 0, 1] = -numpy.cos(E0)*numpy.sin(E2)+numpy.sin(E0)*numpy.sin(E1)*numpy.cos(E2)
    R_old[ 0, 2] =  numpy.sin(E0)*numpy.sin(E2)+numpy.cos(E0)*numpy.sin(E1)*numpy.cos(E2)
    R_old[ 1, 0] = numpy.cos(E1)*numpy.sin(E2)
    R_old[ 1, 1] = numpy.cos(E0)*numpy.cos(E2)+numpy.sin(E0)*numpy.sin(E1)*numpy.sin(E2)
    R_old[ 1, 2] = -numpy.sin(E0)*numpy.cos(E2)+numpy.cos(E0)*numpy.sin(E1)*numpy.sin(E2)
    R_old[ 2, 0] = -numpy.sin(E1)
    R_old[ 2, 1] = numpy.sin(E0)*numpy.cos(E1)
    R_old[ 2, 2] = numpy.cos(E0)*numpy.cos(E1)
    return R_old

def _test_rotation():
    for i in range(10):
        v = numpy.random.rand(3)
        E0 = numpy.random.rand()
        v1 = rotation()
    

def crossproduct(a, b):
    c = numpy.array([a[1]*b[2] - a[2]*b[1],
                     a[2]*b[0] - a[0]*b[2],
                     a[0]*b[1] - a[1]*b[0]])
    return c

def dotproduct(v1, v2):
    return numpy.sum((a*b) for a, b in zip(v1, v2))

def length(v):
    return numpy.sqrt(dotproduct(v, v))

def angle(v1, v2):
    return numpy.arccos(dotproduct(v1, v2) / (length(v1) * length(v2)))

