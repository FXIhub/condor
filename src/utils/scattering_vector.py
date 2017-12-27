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
import numpy, copy

import logging
logger = logging.getLogger(__name__)

from .log import log_and_raise_error,log_warning,log_info,log_debug
#import rotation
import condor.utils.linalg


def q_from_p(p, wavelength):
    r"""
    Return scattering vector from pixel position vector and photon wavelength.

    Args:
      :p (array): Pixel position vector :math:`\vec{p}=(p_x,p_y,p_z)` with respect to the interaction point in unit meter

      :wavelength (float): Photon wavelength in unit meter
    """
    p0 = p / condor.utils.linalg.length(p)
    R_Ewald = 2*numpy.pi / wavelength
    k0 = R_Ewald * numpy.array([0.,0.,1.])
    k1 = R_Ewald * p0
    q = k0 - k1
    return q


def generate_qmap(X,Y,pixel_size,detector_distance,wavelength,extrinsic_rotation=None, order="xyz"):
    r"""
    Generate scattering vector map from experimental parameters

    Args:
      :X (array): :math:`x`-coordinates of pixels in unit meter

      :Y (array): :math:`y`-coordinates of pixels in unit meter

      :pixel_size (float): Pixel size (i.e. edge length) in unit meter

      :detector_distance (float): Distance from interaction point to detector plane

      :wavelength (float): Photon wavelength in unit meter

    Kwargs:
      :extrinsic_rotation (:class:`condor.utils.rotation.Rotation`): Extrinsic rotation of the sample. If ``None`` no rotation is applied (default ``None``)

      :order (str): Order of scattering vector coordinates in the output array. Choose either ``'xyz'`` or ``'zyx'`` (default ``'xyz'``)    
    """
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

def generate_qmap_3d(qn, qmax, extrinsic_rotation=None, order='xyz'):
    q = numpy.linspace(-qmax, qmax, qn)
    Qz, Qy, Qx = numpy.meshgrid(q, q, q, indexing='ij')
    qmap = numpy.zeros(shape=(qn, qn, qn, 3), dtype='float')
    if order == 'xyz':
        qmap[:, :, :, 0] = Qx[:, :, :]
        qmap[:, :, :, 1] = Qy[:, :, :]
        qmap[:, :, :, 2] = Qz[:, :, :]
    elif order == 'zyx':
        qmap[:, :, :, 2] = Qx[:, :, :]
        qmap[:, :, :, 1] = Qy[:, :, :]
        qmap[:, :, :, 0] = Qz[:, :, :]
    else:
        log_and_raise_error(logger, "order=\'%s\' is not a recognised argument for this function." % str(order))
        return
    if extrinsic_rotation is not None:
        log_debug(logger, "Applying qmap rotation.")
        intrinsic_rotation = copy.deepcopy(extrinsic_rotation)
        intrinsic_rotation.invert()
        qmap = intrinsic_rotation.rotate_vectors(qmap.ravel(), order=order).reshape(qmap.shape)
    return qmap

def generate_rpix_3d(qn, qmax, wavelength, detector_distance, pixel_size):
    R_Ewald = 2*numpy.pi/wavelength
    qmap = generate_qmap_3d(qn, qmax)
    q = (qmap**2).sum(axis=3)
    rpix = detector_distance * numpy.tan(2.*numpy.arcsin(qmap/(2.*R_Ewald)))/pixel_size
    return rpix

# Convenience function
def generate_absqmap(X,Y,pixel_size,detector_distance,wavelength):
    """
    Generate absolute scattering vector map from experimental parameters

    Args:
      :X (array): See :func:`condor.utils.scattering_vector.generate_qmap`

      :Y (array): See :func:`condor.utils.scattering_vector.generate_qmap`

      :pixel_size (float): See :func:`condor.utils.scattering_vector.generate_qmap`

      :detector_distance (float): See :func:`condor.utils.scattering_vector.generate_qmap`

      :wavelength (float): See :func:`condor.utils.scattering_vector.generate_qmap`
    """
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

