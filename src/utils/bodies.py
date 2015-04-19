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

import numpy, sys, numpy, types, pickle, time, math
import icosahedron, linalg
 
import logging
logger = logging.getLogger("Condor")
from log import log


def make_icosahedron_map(N,nRmax,euler1=0.,euler2=0.,euler3=0.):
    logger.debug("Building icosahedral geometry")
    logger.debug("Grid: %i x %i x %i (%i voxels)" % (N,N,N,N**3))
    t0 = time.time()
    icomap = icosahedron.icosahedron(N,nRmax,(euler1,euler2,euler3))
    t1 = time.time()
    logger.debug("Built map within %f seconds." % (t1-t0))
    return icomap

def make_icosahedron_map_old(N,nRmax,euler1=0.,euler2=0.,euler3=0.):
    na = nRmax/numpy.sqrt(10.0+2*numpy.sqrt(5))*4.
    nRmin = numpy.sqrt(3)/12*(3.0+numpy.sqrt(5))*na # radius at faces
    logger.debug("Building icosahedral geometry")
    n_list = get_icosahedron_normal_vectors(euler1,euler2,euler3)
    X,Y,Z = 1.0*numpy.mgrid[0:N,0:N,0:N]
    X = X - (N-1)/2.
    Y = Y - (N-1)/2.
    Z = Z - (N-1)/2.
    logger.debug("Grid: %i x %i x %i (%i voxels)" % (N,N,N,N**3))
    icomap = numpy.zeros((len(n_list),N,N,N))
    # calculate distance of all voxels to all faces (negative inside, positive outside icosahedron)
    for i in range(len(n_list)):
        icomap[i,:,:,:] = (X*n_list[i][2]+Y*n_list[i][1]+Z*n_list[i][0])+nRmin
    s = 1.
    M = icomap.copy()
    temp = abs(M)<0.5*s
    icomap[temp] = 0.5+icomap[temp]/s
    icomap[M<(-0.5)*s] = 0
    icomap[M>0.5*s] = 1
    icomap = icomap.min(0)
    return icomap

def get_icosahedron_normal_vectors(euler_1=0.,euler_2=0.,euler_3=0.):
    # construct normal vectors of faces
    phi = (1+numpy.sqrt(5))/2.0
    ri = phi**2/2./numpy.sqrt(3.)
    # normal vectors for every vertice
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
    X = [x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12]
    # angle between normals
    an = round(numpy.dot(x5,x1))

    def cont_element(el,l):
        for i in range(0,len(l)):
            if (el == l[i]).all():
                return True
        return False

    def angles_match(y1,y2,y3):
        if round(numpy.dot(y1,y2)) == an and round(numpy.dot(y2,y3)) == an and round(numpy.dot(y3,y1)) == an:
            return True
        else:
            return False

    n_list = []
    for i in range(0,len(X)):
        for j in range(0,len(X)):
            for k in range(0,len(X)):
                n = (X[i]+X[j]+X[k])/6./ri
                if angles_match(X[i],X[j],X[k]) and not cont_element(n,n_list):
                    n_list.append(n)

                
    if euler_1 != 0. or euler_2 != 0. or euler_3 != 0.:
        for i in range(0,len(n_list)):
            n_list[i] = linalg.rotation(n_list[i],euler_1,euler_2,euler_3)


    return n_list

def make_spheroid_map(N,nA,nB,euler0=0.,euler1=0.,euler2=0.):
    X,Y,Z = 1.0*numpy.mgrid[0:N,0:N,0:N]
    X = X-(N-1)/2.
    Y = Y-(N-1)/2.
    Z = Z-(N-1)/2.
    R_sq = X**2+Y**2+Z**2
    e_c = linalg.rotation(numpy.array([0.0,1.0,0.0]),euler0,euler1,euler2)
    d_sq_c = ((Z*e_c[0])+(Y*e_c[1])+(X*e_c[2]))**2
    r_sq_c = abs( R_sq * (1 - (d_sq_c/(R_sq+numpy.finfo("float32").eps))))
    spheroidmap = r_sq_c/nA**2+d_sq_c/nB**2
    spheroidmap[spheroidmap<=1] = 1
    spheroidmap[spheroidmap>1] = 0
    return spheroidmap

def make_sphere_map(N,nR):
    X,Y,Z = 1.0*numpy.mgrid[0:N,0:N,0:N]
    X = X-(N-1)/2.
    Y = Y-(N-1)/2.
    Z = Z-(N-1)/2.
    R = numpy.sqrt(X**2+Y**2+Z**2)
    spheremap = numpy.zeros(shape=R.shape,dtype="float64")
    spheremap[R<=nR] = 1
    spheremap[abs(nR-R)<0.5] = 0.5+0.5*(nR-R[abs(nR-R)<0.5])
    return spheremap

# A1 is added to A2
def array_to_array(A1,A2,p0=None,origin="corner",mode="sum",fill_value=0.,factor=1.):
    N1 = numpy.array(A1.shape)
    N2 = numpy.array(A2.shape)
    d = len(N1)
    if d > 3:
        logger.error("Cannot handle more than 3 dimensional data.")
        return
    if p0 == None:
        p1 = numpy.zeros(d)
    else:
        p1 = p0
    if origin == "corner":
        p = p1
    elif origin == "middle":
        p = p1+(N2-1)/2.
    else:
        p = p1+origin
    p_min = numpy.int16((p-N1/2).round())
    p_max = p_min + N1
    #print p_min,p_max
    N2_new = N2.copy()
    origin_offset = numpy.zeros(d)
    for di in range(d):
        if p_min[di] < 0:
            offset = -p_min[di]
            N2_new[di] += offset
            p_min[di] += offset
            p_max[di] += offset
        if p_max[di] > N2[di]:
            N2_new[di] = p_max[di]
    A2_new = A2.copy()
    #A2_new = numpy.zeros(shape=tuple(N2_new),dtype=A2.dtype) + fill_value
    #print p_min,p_max
    if mode == "sum": f = lambda a,b: a+b
    elif mode == "replace": f = lambda a,b: b
    elif mode == "max": f = lambda a,b: (a>=b)*a+(b>a)*b
    elif mode == "min": f = lambda a,b: (a<=b)*a+(b<a)*b
    elif mode == "factor": f = lambda a,b: a*(1-b)+b*factor
    else: logger.error("%s is not a valid mode." % mode)
    if d == 1:
        A2_new[p_min[0]:p_max[0]] = f(A2_new[p_min[0]:p_max[0]],A1[:])
    elif d == 2:
        A2_new[p_min[0]:p_max[0],p_min[1]:p_max[1]] = f(A2_new[p_min[0]:p_max[0],p_min[1]:p_max[1]],A1[:,:])
    elif d == 3:
        A2_new[p_min[0]:p_max[0],p_min[1]:p_max[1],p_min[2]:p_max[2]] = f(A2_new[p_min[0]:p_max[0],p_min[1]:p_max[1],p_min[2]:p_max[2]],A1[:,:,:])
    return A2_new
        
