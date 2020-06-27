# -----------------------------------------------------------------------------------------------------
# CONDOR
# Simulator for diffractive single-particle imaging experiments with X-ray lasers
# http://xfel.icm.uu.se/condor/
# -----------------------------------------------------------------------------------------------------
# Copyright 2016 Max Hantke, Filipe R.N.C. Maia, Tomas Ekeberg
# Condor is distributed under the terms of the BSD 2-Clause License
# -----------------------------------------------------------------------------------------------------
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# 
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# -----------------------------------------------------------------------------------------------------
# General note:
# All variables are in SI units by default. Exceptions explicit by variable name.
# -----------------------------------------------------------------------------------------------------

from __future__ import print_function, absolute_import # Compatibility with python 2 and 3
import numpy, sys, numpy, types, pickle, time, math
import condor.utils.icosahedron as icosahedron
import condor.utils.linalg as linalg
 
import logging
logger = logging.getLogger(__name__)
from .log import log_and_raise_error,log_warning,log_info,log_debug


def make_sphere_map(N,nR):
    """
    Generate a 3D map of a sphere particle on a regular grid (values between 0 and 1)

    The result is quite rough (i.e. linear interpolation)

    Args:
      :N (int): Edge length of the grid in unit pixels

      :nR (float): Radius in unit pixels

    .. note:: This function was written for testing purposes and generates a map with rough edges. Use :class:`condor.particle.particle_sphere.ParticleSphere` for more accurate uniform sphere diffraction simulations.
    """
    X,Y,Z = 1.0*numpy.mgrid[0:N,0:N,0:N]
    X = X-(N-1)/2.
    Y = Y-(N-1)/2.
    Z = Z-(N-1)/2.
    R = numpy.sqrt(X**2+Y**2+Z**2)
    spheremap = numpy.zeros(shape=R.shape,dtype="float64")
    spheremap[R<=nR] = 1
    # Linear interpolation at the transition
    spheremap[abs(nR-R)<0.5] = 0.5+0.5*(nR-R[abs(nR-R)<0.5])
    return spheremap


def make_spheroid_map(N, nA, nC, rotation=None):
    """
    Generate a 3D binary map of a spheroid particle on a regular grid

    The result is very rough (i.e. nearest-neighbor interpolation)

    Args:
      :N (int): Edge length of the grid in unit pixels

      :nA (float): Radius perpendicular to the rotation axis of the ellipsoid in unit pixels

      :nC (float): Radius along the rotation axis of the ellipsoid in unit pixels

    Kwargs:
    
      :rotation (:class:`condor.utils.rotation.Rotation`): Rotation instance for extrinsic rotation of the icosahedron. 

    .. note:: This function was written for testing purposes and generates a map with rough edges. Use :class:`condor.particle.particle_spheroid.ParticleSpheroid` for more accurate uniform spheroid diffraction simulations.
    """
    X,Y,Z = 1.0*numpy.mgrid[0:N,0:N,0:N]
    X = X-(N-1)/2.
    Y = Y-(N-1)/2.
    Z = Z-(N-1)/2.
    R_sq = X**2+Y**2+Z**2
    e_c = numpy.array([0.0,1.0,0.0])
    if rotation is not None:
        e_c = rotation.rotate_vector(e_c)
    d_sq_c = ((X*e_c[0])+(Y*e_c[1])+(Z*e_c[2]))**2
    r_sq_c = abs( R_sq * (1 - (d_sq_c/(R_sq+numpy.finfo("float32").eps))))
    spheroidmap = r_sq_c/float(nA)**2+d_sq_c/float(nC)**2
    spheroidmap[spheroidmap<=1] = 1
    spheroidmap[spheroidmap>1] = 0
    return spheroidmap

def make_icosahedron_map(N,nRmax,extrinsic_rotation=None):
    """
    Generate map of a uniform icosahedron (density = 1) on a regular grid

    Orientation: The cartesian grid axis all lie parallel to 2-fold symmetry axes of the icosahedron.

    Args:
      :N (int): Edge length of the grid in unit pixels

      :nRmax (float): Outer radius of the icosahedron in unit pixels

    Kwargs:
      :rotation (:class:`condor.utils.rotation.Rotation`): Rotation instance for extrinsic rotation of the icosahedron. 
    """
    log_debug(logger, "Building icosahedral geometry")
    log_debug(logger, "Grid: %i x %i x %i (%i voxels)" % (N,N,N,N**3))
    t0 = time.time()
    if extrinsic_rotation is not None:
        q = extrinsic_rotation.get_as_quaternion()
        icomap = icosahedron.icosahedron(N,nRmax,q)
    else:
        icomap = icosahedron.icosahedron(N,nRmax)
    t1 = time.time()
    log_debug(logger, "Built map within %f seconds." % (t1-t0))
    return icomap

def make_icosahedron_map_slow(N,nRmax,extrinsic_rotation=None):
    """
    Generate map of a uniform icosahedron (density = 1) on a regular grid (*slow python implementation*)

    Orientation: The cartesian grid axis all lie parallel to 2-fold symmetry axes of the icosahedron.

    Args:
      :N (int): Edge length of the grid in unit pixels

      :nRmax (float): Outer radius of the icosahedron in unit pixels

    Kwargs:
    
      :rotation (:class:`condor.utils.rotation.Rotation`): Rotation instance for extrinsic rotation of the icosahedron. 
    """
    na = nRmax/numpy.sqrt(10.0+2*numpy.sqrt(5))*4.
    nRmin = numpy.sqrt(3)/12*(3.0+numpy.sqrt(5))*na # radius at faces
    log_debug(logger, "Building icosahedral geometry")
    n_list = get_icosahedron_normal_vectors()
    # Rotate
    if extrinsic_rotation is not None:
        n_list = extrinsic_rotation.rotate_vectors(numpy.array(n_list))
    X,Y,Z = 1.0*numpy.mgrid[0:N,0:N,0:N]
    X = X - (N-1)/2.
    Y = Y - (N-1)/2.
    Z = Z - (N-1)/2.
    log_debug(logger, "Grid: %i x %i x %i (%i voxels)" % (N,N,N,N**3))
    icomap = numpy.zeros((len(n_list),N,N,N))
    # calculate distance of all voxels to all faces (negative inside, positive outside icosahedron)
    for i in range(len(n_list)):
        icomap[i,:,:,:] = (X*n_list[i][0]+Y*n_list[i][1]+Z*n_list[i][2])+nRmin
    s = 1.
    M = icomap.copy()
    temp = abs(M)<0.5*s
    icomap[temp] = 0.5+icomap[temp]/s
    icomap[M<(-0.5)*s] = 0
    icomap[M>0.5*s] = 1
    icomap = icomap.min(0)
    return icomap

def get_icosahedron_vertices():
    """
    Return array of vertices vectors of a regular icosahedron
    """
    # Weisstein, Eric W. "Icosahedral Group." From MathWorld--A Wolfram Web Resource. http://mathworld.wolfram.com/IcosahedralGroup.html
    phi = (1+numpy.sqrt(5))/2.0
    x1 = numpy.array([0.0,1.0,phi])
    x2 = numpy.array([0.0,1.0,-phi])
    x3 = numpy.array([0.0,-1.0,phi])
    x4 = numpy.array([0.0,-1.0,-phi])
    x5 = numpy.array([1.0,phi,0.0])
    x6 = numpy.array([1.0,-phi,0.0])
    x7 = numpy.array([-1.0,phi,0.0])
    x8 = numpy.array([-1.0,-phi,0.0])
    x9 = numpy.array([phi,0.0,1.0])
    x10 = numpy.array([-phi,0.0,1.0])
    x11 = numpy.array([phi,0.0,-1.0])
    x12 = numpy.array([-phi,0.0,-1.0])
    return numpy.array([x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12])
    
def get_icosahedron_normal_vectors():
    """
    Return array of normal vectors at the faces of a regular icosahedron
    """

    X = get_icosahedron_vertices()
    # Inner radius of icosahedron (distance from center to each face)
    # Weisstein, Eric W. "Icosahedral Group." From MathWorld--A Wolfram Web Resource. http://mathworld.wolfram.com/IcosahedralGroup.html
    phi = (1+numpy.sqrt(5))/2.0
    ri = phi**2/2./numpy.sqrt(3.)
    # Convenience functions
    # Check whether vector is in list
    def cont_element(el,l):
        for i in range(0,len(l)):
            if (el == l[i]).all():
                return True
        return False
    # Angle between neighboring vertices
    an = round(numpy.dot(X[4],X[0]),1)
    # Check whether vertices are neighbors
    def neighbors(v1,v2,v3):
        if round(numpy.dot(v1,v2),1) == an and round(numpy.dot(v2,v3),1) == an and round(numpy.dot(v3,v1),1) == an:
            return True
        else:
            return False
    n_list = []
    # Loop through all 3-vector permutations and build normal vector from the neighboring vertices
    for i in range(0,len(X)):
        for j in range(0,len(X)):
            for k in range(0,len(X)):
                n = (X[i]+X[j]+X[k])/6./ri
                if neighbors(X[i],X[j],X[k]) and not cont_element(n,n_list):
                    n_list.append(n)
    return n_list

