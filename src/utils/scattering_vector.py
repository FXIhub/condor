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

import numpy, copy

import logging
logger = logging.getLogger(__name__)

from log import log_and_raise_error,log_warning,log_info,log_debug
import rotation
import linalg


def q_from_p(p, wavelength):
    p0 = p / linalg.length(p)
    R_Ewald = 2*numpy.pi / wavelength
    k0 = R_Ewald * numpy.array([0.,0.,1.])
    k1 = R_Ewald * p0
    q = k0 - k1
    return q

def generate_qmap(X,Y,pixel_size,detector_distance,wavelength,extrinsic_rotation=None, order="xyz"):
    log_debug(logger, "Allocating qmap (%i,%i,%i)" % (X.shape[0],X.shape[1],3))
    R_Ewald = 2*numpy.pi/wavelength
    p_x = X*pixel_size
    p_y = Y*pixel_size
    p_z = numpy.ones_like(X)*detector_distance
    l = numpy.sqrt(p_x**2+p_y**2+p_z**2)
    r_x = p_x/l
    r_y = p_y/l
    r_z = p_z/l - 1.
    qmap = numpy.zeros(shape=(X.shape[0],X.shape[1],3))
    if order == "xyz":
        qmap[:,:,0] = r_x * R_Ewald
        qmap[:,:,1] = r_y * R_Ewald
        qmap[:,:,2] = r_z * R_Ewald
    elif order == "zyx":
        qmap[:,:,0] = r_z * R_Ewald
        qmap[:,:,1] = r_y * R_Ewald
        qmap[:,:,2] = r_x * R_Ewald
    else:
        log_and_raise_error(logger, "Indexing with order=%s is invalid." % order)
    if extrinsic_rotation is not None:
        log_debug(logger, "Applying qmap rotation.")
        intrinsic_rotation = copy.deepcopy(extrinsic_rotation)
        intrinsic_rotation.invert()
        qmap = intrinsic_rotation.rotate_vectors(qmap.ravel(), order=order).reshape(qmap.shape)
    return qmap

def generate_absqmap(X,Y,pixel_size,detector_distance,wavelength):
    qmap = generate_qmap(X,Y,pixel_size,detector_distance,wavelength)
    qmap = numpy.sqrt(qmap[:,:,0]**2+qmap[:,:,1]**2+qmap[:,:,2]**2)
    return qmap

#def generate_qmap_ori(X,Y,pixel_size,detector_distance,wavelength):
#    phi = numpy.arctan2(pixel_size*numpy.sqrt(X**2+Y**2),detector_distance)
#    R_Ewald = 2*numpy.pi/wavelength
#    qx = R_Ewald*2*numpy.sin(numpy.arctan2(pixel_size*X,detector_distance)/2.)
#    qy = R_Ewald*2*numpy.sin(numpy.arctan2(pixel_size*Y,detector_distance)/2.)
#    qz = R_Ewald*(1-numpy.cos(phi))
#    qmap = numpy.zeros(shape=(X.shape[0],Y.shape[1],3))
#    qmap[:,:,0] = qx[:,:]
#    qmap[:,:,1] = qy[:,:]
#    qmap[:,:,2] = qz[:,:]
#    return qmap

