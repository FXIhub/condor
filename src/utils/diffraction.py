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
