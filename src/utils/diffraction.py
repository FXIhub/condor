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
    return wavelength / 2. / numpy.arcsin( numpy.tan( pixel_center_distance / detector_distance ) / 2.)

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
