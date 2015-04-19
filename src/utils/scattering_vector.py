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
 
import logging
logger = logging.getLogger("Condor")
from log import log

def generate_qmap(X,Y,pixel_size,detector_distance,wavelength,euler_angle_0=0.,euler_angle_1=0.,euler_angle_2=0.):
    log(logger.debug,"Allocating qmap.")
    R_Ewald = 2*numpy.pi/wavelength
    qx = R_Ewald*2*numpy.sin(numpy.arctan2(pixel_size*X,detector_distance)/2.)
    qy = R_Ewald*2*numpy.sin(numpy.arctan2(pixel_size*Y,detector_distance)/2.)
    phi = numpy.arctan2(pixel_size*numpy.sqrt(X**2+Y**2),detector_distance)
    qz = R_Ewald*(1-numpy.cos(phi))
    qmap = numpy.zeros(shape=(X.shape[0],Y.shape[1],3))
    qmap[:,:,0] = qz[:,:]
    qmap[:,:,1] = qy[:,:]
    qmap[:,:,2] = qx[:,:]
    if euler_angle_0 != 0. or euler_angle_1 != 0. or euler_angle_2 != 0.:
        log(logger.debug,"Applying qmap rotation with angles %f, %f, %f." % (euler_angle_0,euler_angle_1,euler_angle_2))
        # Old and slow
        #for iy in numpy.arange(0,qmap.shape[0]):
        #    for ix in numpy.arange(0,qmap.shape[1]):
        #        qmap[iy,ix,:] = rotation(qmap[iy,ix,:],euler_angle_0,euler_angle_1,euler_angle_2)
        cos = numpy.cos
        sin = numpy.sin
        E0 = euler_angle_0
        E1 = euler_angle_1
        E2 = euler_angle_2
        M = numpy.array([[cos(E1)*cos(E2),
                          -cos(E0)*sin(E2)+sin(E0)*sin(E1)*cos(E2),
                          sin(E0)*sin(E2)+cos(E0)*sin(E1)*cos(E2)],
                         [cos(E1)*sin(E2),
                          cos(E0)*cos(E2)+sin(E0)*sin(E1)*sin(E2),
                          -sin(E0)*cos(E2)+cos(E0)*sin(E1)*sin(E2)],
                         [-sin(E1),
                          sin(E0)*cos(E1),
                          cos(E0)*cos(E1)]])
        Y,X = numpy.mgrid[:qmap.shape[0],:qmap.shape[1]]
        Y = Y.flatten()
        X = X.flatten()
        s = qmap.shape
        qmap = numpy.array([numpy.dot(M,qmap[iy,ix,:]) for ix,iy in zip(X,Y)])
        qmap = qmap.reshape(s)
    return qmap

def generate_absqmap(X,Y,pixel_size,detector_distance,wavelength):
    qmap = generate_qmap(X,Y,pixel_size,detector_distance,wavelength)
    qmap = numpy.sqrt(qmap[:,:,0]**2+qmap[:,:,1]**2+qmap[:,:,2]**2)
    return qmap

def generate_qmap_ori(X,Y,pixel_size,detector_distance,wavelength):
    phi = numpy.arctan2(pixel_size*numpy.sqrt(X**2+Y**2),detector_distance)
    R_Ewald = 2*numpy.pi/wavelength
    qx = R_Ewald*2*numpy.sin(numpy.arctan2(pixel_size*X,detector_distance)/2.)
    qy = R_Ewald*2*numpy.sin(numpy.arctan2(pixel_size*Y,detector_distance)/2.)
    qz = R_Ewald*(1-numpy.cos(phi))
    qmap = numpy.zeros(shape=(X.shape[0],Y.shape[1],3))
    qmap[:,:,0] = qz[:,:]
    qmap[:,:,1] = qy[:,:]
    qmap[:,:,2] = qx[:,:]
    return qmap

x_to_q = lambda x,pixel_size,detector_distance,wavelength: 4*numpy.pi/wavelength*numpy.sin(numpy.arctan(x*pixel_size/detector_distance)/2.)
