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
 
import logging
logger = logging.getLogger("Condor")
from log import log

def random_euler_angles():
    """
    Generates a triplet (phi, theta, psi) of random Euler angles.
    """
    r1,r2,r3 = numpy.random.random(3)
    q1 = numpy.sqrt(1.0-r1)*numpy.sin(2.0*numpy.pi*r2)
    q2 = numpy.sqrt(1.0-r1)*numpy.cos(2.0*numpy.pi*r2)
    q3 = numpy.sqrt(r1)*numpy.sin(2.0*numpy.pi*r3)
    q4 = numpy.sqrt(r1)*numpy.cos(2.0*numpy.pi*r3)
    e1 = math.atan2(2.0*(q1*q2+q3*q4), 1.0-2.0*(q2**2+q3**2))
    e2 = math.asin(2.0*(q1*q3-q4*q2))
    e3 = math.atan2(2.0*(q1*q4+q2*q3), 1.0-2.0*(q3**2+q4**2))
    return (e1,e2,e3)

def crystallographic_resolution(wavelength, pixel_center_distance, detector_distance):
    """
    Returns crystallographic resolution (full-period resolution at the closest edge)
    """
    return wavelength / numpy.sin( numpy.arctan( pixel_center_distance / detector_distance ) )
    
def nyquist_pixel_size(wavelength, detector_distance, particle_size):
    """
    Returns size of one Nyquist pixel on the detector in meter.
    """
    return detector_distance * wavelength / particle_size
