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

import numpy

import scipy.constants as constants

import logging
logger = logging.getLogger(__name__)

from log import log_and_raise_error,log_warning,log_info,log_debug

class Photon:
    """
    Class for X-ray photon
    """
    def __init__(self, wavelength=None, energy=None, energy_eV=None):
        """
        Initialisation of a Photon instance

        The photon can be initialised either by passing the wavelength, the photon energy in unit Joule or the photon energy in unit eV

        Kwargs:
           :wavelength(float): Photon wavelength [m]
           :energy(float): Photon energy [J]
           :energy_eV(float): Photon energy in electron volts [eV]
        """
        if (wavelength is not None and energy is not None) or (wavelength is not None and energy_eV is not None) or (energy is not None and energy_eV is not None):
            log_and_raise_error(logger, "Invalid arguments during initialisation of Photon instance. More than one of the arguments is not None.")
            return
        if wavelength is not None:
            self.set_wavelength(wavelength)
        elif energy is not None:
            self.set_energy(energy,"J")
        elif energy_eV is not None:
            self.set_energy(energy_eV,"eV")
        else:
            log_and_raise_error(logger, "Photon could not be initialized. It needs to be initialized with either the wavelength, the photon energy in unit Joule or the photon energy in unit eV.")
            
    def get_energy(self,unit="J"):
        if unit == "J":
            return self._energy
        elif unit == "eV":
            return self._energy/constants.e
        else:
            log_and_raise_error(logger, "%s is not a valid energy unit." % unit)

    def set_energy(self,energy,unit="J"):
        if unit == "J":
            self._energy = energy
        elif unit == "eV":
            self._energy = energy*constants.e
        else:
            log_and_raise_error(logger, "%s is not a valid energy unit." % unit)

    def get_wavelength(self):
        return constants.c*constants.h/self._energy

    def set_wavelength(self,wavelength):
        self._energy = constants.c*constants.h/wavelength

        
