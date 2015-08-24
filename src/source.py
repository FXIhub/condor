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

import sys,os
import numpy

import scipy.constants as constants

import logging
logger = logging.getLogger(__name__)

import condor.utils.log
from condor.utils.log import log_and_raise_error,log_warning,log_info,log_debug
from condor.utils.config import load_config
from condor.utils.variation import Variation

def load_source(conf):
    """
    Create new Source instance and load parameters from a Condor configuration file.
    
    Args:
       :conf(str): Condor configuration file
    """
    C = condor.utils.config.load_config({"source": load_config(conf)["source"]}, {"source": load_config(condor.CONDOR_default_conf)["source"]})
    source = Source(**C["source"])
    return source
  
class Source:
    """
    Class for the X-ray source
    """
    def __init__(self, wavelength, focus_diameter, pulse_energy, profile_model=None, pulse_energy_variation=None, pulse_energy_spread=None, pulse_energy_variation_n=None, number_of_shots=1):
        """
        Initialisation of a Source instance

        Args:
            :wavelength(float): X-ray wavelength [m]
            :focus_diameter(float): Focus diameter (characteristic transverse dimension) [m]
            :pulse_energy(float): Pulse energy (or its statistical mean) [J]

        Kwargs:
            :number_of_shots(int): Number of independent shots to be simulated (default = 1)
            :pulse_energy_variation(str): Statistical variation of the pulse energy - either \"normal\", \"uniform\", \"range\" or None (no variation of the pulse energy) (default = None)
            :pulse_energy_spread(float): Statistical spread of the pulse energy [J] (default = None)
            :pulse_energy_variation_n(int): This parameter is only effective if pulse_energy_variation=\"range\". In that case this parameter determines the number of samples within the interval pulse_energy +/- pulse_energy_spread/2 (default = None)
            :profile_model(str): Model for the spatial illumination profile - either \"top_hat\", \"pseudo_lorentzian\", \"gaussian\" or None (infinite extent of the illuminatrion, intensity same as in case of \"top_hat\") (default = None)
        """
        self.photon = Photon(wavelength=wavelength)
        self.pulse_energy_mean = pulse_energy
        self.set_pulse_energy_variation(variation=pulse_energy_variation, spread=pulse_energy_spread, n=pulse_energy_variation_n)
        self.profile = Profile(model=profile_model, focus_diameter=focus_diameter)
        self.number_of_shots = number_of_shots
        log_debug(logger, "Source configured")

    def get_conf(self):
        conf = {}
        conf["source"] = {}
        conf["source"]["wavelength"]               = self.photon.get_wavelength()
        conf["source"]["focus_diameter"]           = self.profile.focus_diameter
        conf["source"]["pulse_energy"]             = self.pulse_energy_mean
        conf["source"]["profile_model"]            = self.profile.get_model()
        pevar = self._pulse_energy_variation.get_conf()
        conf["source"]["pulse_energy_variation"]   = pevar["mode"]
        conf["source"]["pulse_energy_spread"]      = pevar["spread"]
        conf["source"]["pulse_energy_variation_n"] = pevar["n"]
        conf["source"]["number_of_shots"]          = self.number_of_shots
        return conf
        
    def set_pulse_energy_variation(self, variation = None, spread = None, n = None):
        """
        Specify how the pulse energy varies. If no inputs are given the pulse energy is constant

        Kwargs:
            :pulse_energy_variation(str): Statistical variation of the pulse energy - either \"normal\", \"uniform\", \"range\" or None (no variation of the pulse energy) (default = None)
            :pulse_energy_spread(float): Statistical spread of the pulse energy [J] (default = None)
            :pulse_energy_variation_n(int): This parameter is only effective if pulse_energy_variation=\"range\". In that case this parameter determines the number of samples within the interval pulse_energy +/- pulse_energy_spread/2 (default = None)
        """
        self._pulse_energy_variation = Variation(variation,spread,n,number_of_dimensions=1,name="pulse energy")

    def get_intensity(self, position, unit = "ph/m2", pulse_energy = None):
        """
        Retrieve the intensity at a given position in the beam

        Args:
           :position(list): Coordinates of the position where the intensity shall be calculated
           
        Kwargs:
           :unit(str): Intensity unit - either \"ph/m2\", \"J/m2\", \"J/um2\" or \"mJ/um2\" (default = \"ph/m2\")
           :pulse_energy(float): Pulse energy - if None the statistical mean of the pulse energy will be used [J] (default = None)
        """
        # Assuming
        # 1) Radially symmetric profile that is invariant along the beam axis within the sample volume
        # 2) The variation of intensity are on much larger scale than the dimension of the particle size (i.e. flat wavefront)
        r = numpy.sqrt(position[1]**2 + position[2]**2)
        I = (self.profile.get_radial())(r) * (pulse_energy if pulse_energy is not None else self.pulse_energy_mean)
        if unit == "J/m2":
            pass
        elif unit == "ph/m2":
            I /= self.photon.get_energy("J") 
        elif unit == "J/um2":
            I *= 1.E-12
        elif unit == "mJ/um2":
            I *= 1.E-9
        else:
            log_and_raise_error(logger, "%s is not a valid unit." % unit)
            return
        return I

    def get_next(self):
        """
        Iterate and return parameters in a dictionary
        """
        return {"pulse_energy":self._get_next_pulse_energy(),
                "wavelength":self.photon.get_wavelength(),
                "photon_energy":self.photon.get_energy("J"),
                "photon_energy_eV":self.photon.get_energy("eV")}

    def _get_next_pulse_energy(self):
        p = self._pulse_energy_variation.get(self.pulse_energy_mean)
        # Non-random
        if self._pulse_energy_variation._mode in [None,"range"]:
            if p <= 0:
                log_and_raise_error(logger, "Pulse energy smaller-equals zero. Change your configuration.")
            else:
                return p
        # Random
        else:
            if p <= 0.:
                log_warning(logger, "Pulse energy smaller-equals zero. Try again.")
                self._get_next_pulse_energy()
            else:
                return p


    
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

        
class Profile:
    """
    Class for spatial illumination profile
    """
    def __init__(self, model, focus_diameter):
        """
        Initialisation of a Profile instance

        Args:
           :model(str): Model for the spatial illumination profile - either \"top_hat\", \"pseudo_lorentzian\", \"gaussian\" or None (infinite extent of the illuminatrion, intensity same as in case of \"top_hat\")
           :focus_diameter(float): Focus diameter (or characteristic dimension) [m]
        """
        self.set_model(model)
        self.focus_diameter = focus_diameter
        
    def set_model(self,model):
        if model is None or model in ["top_hat","pseudo_lorentzian","gaussian"]:
            self._model = model
        else:
            log_and_raise_error(logger, "Pulse profile model %s is not implemented. Change your configuration and try again.")
            sys.exit(0)

    def get_model(self):
        return self._model

    def get_radial(self):
        if self._model is None:
            # we always hit with full power
            p = lambda r: 1. / (numpy.pi * (self.focus_diameter / 2.)**2)
        elif self._model == "top_hat":
            # focus diameter is diameter of top hat profile
            def p(r):
                if numpy.isscalar(r):
                    return (1.  / (numpy.pi * (self.focus_diameter / 2.)**2)) if r < (self.focus_diameter / 2.) else 0.
                else:
                    return (1.  / (numpy.pi * (self.focus_diameter / 2.)**2)) * (r < (self.focus_diameter / 2.))
        elif self._model == "pseudo_lorentzian":
            # focus diameter is FWHM of lorentzian
            sigma = self.focus_diameter / 2.
            p = lambda r: _pseudo_lorentzian_2dnorm(r, sigma)
        elif self._model == "gaussian":
            # focus diameter is FWHM of gaussian
            sigma = self.focus_diameter / (2.*numpy.sqrt(2.*numpy.log(2.)))
            p = lambda r: _gaussian_2dnorm(r, sigma)
        return p
            
_gaussian = lambda x, sigma: numpy.exp(-x**2/(2*sigma**2))

_gaussian_2dnorm = lambda x, sigma: _gaussian(x, sigma) / ( 2 * numpy.pi * sigma**2 )

_lorentzian = lambda x, sigma: sigma**2 / (x**2 + sigma**2)

_pseudo_lorenzian_A1 = 0.74447313315648778 
_pseudo_lorenzian_A2 = 0.22788162774723308
_pseudo_lorenzian_s1 = 0.73985516665883544
_pseudo_lorenzian_s2 = 2.5588165723260907
_pseudo_lorentzian = lambda x, sigma: _pseudo_lorenzian_A1 * _gaussian(x, _pseudo_lorenzian_s1*sigma) + \
                                      _pseudo_lorenzian_A2 * _gaussian(x, _pseudo_lorenzian_s2*sigma)

_pseudo_lorentzian_2dnorm = lambda x, sigma: _pseudo_lorentzian(x, sigma) / ( 2. * numpy.pi * ( _pseudo_lorenzian_A1 * (_pseudo_lorenzian_s1*sigma)**2 + \
                                                                                                _pseudo_lorenzian_A2 * (_pseudo_lorenzian_s2*sigma)**2 ) )

