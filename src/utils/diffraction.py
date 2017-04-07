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

import numpy, math
 
def crystallographic_resolution(wavelength, pixel_center_distance, detector_distance):
    r"""
    Returns crystallographic resolution :math:`R_f` (full-period resolution) in unit meter

    .. math::

      R_f = \frac{ \lambda }{ 2 \sin\left( \arctan\left( \frac{X}{D} \right) / 2 \right) }

    Args:
      :wavelength (float): Photon wavelength :math:`\lambda` in unit meter

      :pixel_center_distance (float): Distance :math:`X` between beam center and pixel measured orthogonally with respect to the beam axis. The unit is meter
    
      :detector_distance: Distance :math:`D` between interaction point and detector plane in unit meter
    """
    return wavelength / 2. / numpy.sin( numpy.arctan( pixel_center_distance / detector_distance ) / 2.)

def resolution_element(wavelength, pixel_center_distance, detector_distance):
    r"""
    Returns length :math:`R_h` of one resolution element (half-period resolution) in unit meter

    .. math::

      R_h = \frac{ \lambda }{ 4 \, \sin\left( \arctan \left( \frac{X}{D} \right) / 2 \right) }

    Args:
      :wavelength (float): Photon wavelength :math:`\lambda` in unit meter

      :pixel_center_distance (float): Distance :math:`X` between beam center and pixel measured orthogonally with respect to the beam axis. The unit is meter
    
      :detector_distance: Distance :math:`D` between interaction point and detector plane in unit meter
    """
    return (0.5*crystallographic_resolution(wavelength, pixel_center_distance, detector_distance))

def nyquist_pixel_size(wavelength, detector_distance, particle_size):
    r"""
    Returns size :math:`p_N` of one Nyquist pixel on the detector in unit meter

    .. math::
    
      p_N = \frac{ D \lambda }{ d }

    Args:
      :wavelength (float): Photon wavelength :math:`\lambda` in unit meter
    
      :detector_distance (float): Distance :math:`D` between interaction point and detector plane in unit meter

      :particle_size (float): Size or characteristic dimension :math:`d` of the particle in unit meter
    """
    return detector_distance * wavelength / particle_size

def polarization_factor(x, y, detector_distance, polarization="ignore"):
    """
    Returns polarization factor for a given geometry and polarization
    
    Horizontally polarized:

    .. math::

      P = \cos^2\left(\arcsin\left(\frac{x}{\sqrt{x^2+y^2+D^2}}\right)\right)

    Vertically polarized:

    .. math::

      P = \cos^2\left(\arcsin\left(\frac{y}{\sqrt{x^2+y^2+D^2}}\right)\right)

    Unpolarized:
    
      P = \left(1 + \cos^2\left(\frac{\sqrt{x^2+y^2}}{\sqrt{x^2+y^2+D^2}}\right)^2\right)

    Ignore polarization:

    .. math::
    
      P = 1

    Args:
      :x (float): horizontal pixel coordinate :math:`x` with respect to beam center in unit meter

      :y (float): vertical pixel coordinate :math:`y` with respect to beam center in unit meter

      :detector_distance (float): detector distance :math:`D` in unit meter

      :polarization (str): Type of polarization can be either *vertical*, *horizontal*, *unpolarized*, or *ignore* 
    """
    if polarization not in ["ignore", "vertical", "horizontal", "unpolarized"]:
        log_and_raise_error(logger, "polarization=\"%s\" is an invalid argument for this function." % polarization)
        return
    if polarization == "ignore":
        P = 1.
    else:
        r = numpy.sqrt(x**2 + y**2 + detector_distance**2)
        if polarization == "vertical":
            P = numpy.cos( numpy.arcsin(y/r) )**2
        elif polarization == "horizontal":
            P = numpy.cos( numpy.arcsin(x/r) )**2
        elif polarization == "unpolarized":
            P = ( 1. + numpy.cos( numpy.arcsin(numpy.sqrt(x**2+y**2)/r) )**2 )
    return P
