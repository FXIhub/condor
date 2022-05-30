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
import numpy
import scipy.constants as constants
import unittest

import logging
logger = logging.getLogger(__name__)

from .log import log_and_raise_error,log_warning,log_info,log_debug

class Photon:
    """
    Class for X-ray photon

    A Photon instance is initialised with one keyword argument - either with ``wavelength``, ``energy`` or ``energy_eV``.
    
    Kwargs:
      :wavelength (float): Photon wavelength in unit meter

      :energy (float): Photon energy in unit Joule

      :energy_eV (float): Photon energy in unit electron volt

      :frequency (float): Photon frequency in unit Hz
    """
    def __init__(self, wavelength=None, energy=None, energy_eV=None, frequency=None):
        if (wavelength is not None and energy is not None) or (wavelength is not None and energy_eV is not None) or (energy is not None and energy_eV is not None):
            log_and_raise_error(logger, "Invalid arguments during initialisation of Photon instance. More than one of the arguments is not None.")
            return
        if wavelength is not None:
            self.set_wavelength(wavelength)
        elif energy is not None:
            self.set_energy(energy)
        elif energy_eV is not None:
            self.set_energy_eV(energy_eV)
        elif frequency is not None:
            self.set_frequency(frequency)
        else:
            log_and_raise_error(logger,
                                "Photon could not be initialised."
                                "It has to be initialized with exactly one of the following keyword arguments"
                                "\t- wavelegth:\t photon wavelength in meters"
                                "\t- energy:\t photon energy in Joules"
                                "\t- energy_eV:\t photon energy in electron volts"
                                "\t- frequency:\t photon frequency in Hertz")
                            
    def set_energy(self, energy):
        """
        Set photon energy in unit Joule

        Args:
          :energy (float): Photon energy in unit Joule
        """
        self._energy = energy

    def set_energy_eV(self, energy_eV):
        """
        Set photon energy in unit electron volt

        Args:
          :energy (float): Photon energy in unit electron volt
        """
        self._energy = energy_eV*constants.e


    def set_wavelength(self, wavelength):
        """
        Set photon wavelength

        Args:
          :wavelength (float): Photon wavelength in unit meter
        """
        self._energy = constants.c*constants.h/wavelength

    def set_frequency(self, frequency):
        """
        Set photon frequency

        Args:
          :frequency (float): Photon frequency in unit Hertz
        """
        self._energy = constants.h*frequency
                                
            
    def get_energy(self):
        """
        Return photon energy in unit Joule
        """
        return self._energy

    def get_energy_eV(self):
        """
        Return photon energy in unit electron volt
        """
        return self._energy/constants.e


    def get_wavelength(self):
        """
        Return wavelength in unit meter
        """
        return constants.c*constants.h/self._energy

    def get_frequency(self):
        """
        Return frequency in unit Hertz
        """
        return self._energy/constants.h

                
        
        
