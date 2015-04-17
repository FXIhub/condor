# ----------------------------------------------------------------------------------------------------- 
# CONDOR 
# Simulator for diffractive single-particle imaging experiments with X-ray lasers
# http://xfel.icm.uu.se/condor/
# ----------------------------------------------------------------------------------------------------- 
# Copyright 2014 Max Hantke, Filipe R.N.C. Maia, Tomas Ekeberg
# Condor is distributed under the terms of the GNU General Public License
# ----------------------------------------------------------------------------------------------------- 
# General note:
#  All variables are in SI units by default. Exceptions explicit by variable name.
# ----------------------------------------------------------------------------------------------------- 

import os
import numpy
import config

import scipy.constants as constants

import logging
logger = logging.getLogger("Condor")
import utils.log
from utils.log import log 

import config
from utils.variation import Variation



class Source:
    """
    A subclass of the input object.
    Defines properties of the FEL x-ray pulse.

    """
    def __init__(self,**kwargs):       

        # Check for valid set of keyword arguments
        req_keys = ["wavelength", "focus_diameter", "pulse_energy"]
        opt_keys = ["pulse_energy_variation", "pulse_energy_spread", "pulse_energy_variation_n", "profile_model"]
        miss_keys,ill_keys = config.check_input(kwargs.keys(),req_keys,opt_keys)
        if len(miss_keys) > 0: 
            for k in miss_keys:
                log(logger.error,"Cannot initialize Source instance. %s is a necessary keyword." % k)
            exit(1)
        if len(ill_keys) > 0:
            for k in ill_keys:
                log(logger.error,"Cannot initialize Source instance. %s is an illegal keyword." % k)
            exit(1)

        # Start initialisation
        self.photon = Photon(wavelength=kwargs["wavelength"])
        if "pulse_energy" in kwargs:
            # Maintain depreciated keyword
            self.pulse_energy_mean = kwargs["pulse_energy"]
        else:
            self.pulse_energy_mean = kwargs["pulse_energy_mean"]
        self.set_pulse_energy_variation(variation=kwargs.get("pulse_energy_variation",None),
                                        spread=kwargs.get("pulse_energy_spread",None),
                                        n=kwargs.get("pulse_energy_variation_n",None))
        self.profile = Profile(model=kwargs.get("profile_model",None),focus_diameter=kwargs["focus_diameter"])
        log(logger.debug,"Source configured")

    def set_pulse_energy_variation(self,variation=None, spread=None, n=None):
        self._pulse_energy_variation = Variation(variation,spread,n,number_of_dimensions=1,name="pulse energy")

    def get_intensity(self,position,unit="ph/m2"):
        # Assuming
        # 1) Radially symmetric profile that is invariant along the beam axis within the sample volume
        # 2) The variation of intensity are on much larger scale than the dimension of the particle size (i.e. flat wavefront)
        r = numpy.sqrt(position[1]**2 + position[2]**2)
        I = (self.profile.get_radial())(r) * self.pulse_energy
        if unit == "J/m2":
            pass
        elif unit == "ph/m2":
            I /= self.photon.get_energy("J") 
        elif unit == "J/um2":
            I *= 1.E-12
        elif unit == "mJ/um2":
            I *= 1.E-9
        else:
            log(logger.error,"%s is not a valid unit." % unit)
            return
        return I

    #def get_area(self):
    #    return numpy.pi*(self.focus_diameter/2.0)**2
                
    def get_next(self):
        self._next_pulse_energy()
        return {"pulse_energy":self.pulse_energy,
                "wavelength":self.photon.get_wavelength()}

    def _next_pulse_energy(self):
        p = self._pulse_energy_variation.get(self.pulse_energy_mean)
        # Non-random
        if self._pulse_energy_variation._mode in [None,"range"]:
            if p <= 0:
                log(logger.error,"Pulse energy smaller-equals zero. Change your configuration.")
            else:
                self.pulse_energy = p
        # Random
        else:
            if p <= 0.:
                log(logger.warning,"Pulse energy smaller-equals zero. Try again.")
                self._next_pulse_energy()
            else:
                self.pulse_energy = p


    
class Photon:
    def __init__(self,**kwarg):
        if "wavelength" in kwarg.keys(): self.set_wavelength(kwarg["wavelength"])
        elif "energy" in kwarg.keys(): self.set_energy(kwarg["energy"],"J")
        elif "energy_eV" in kwarg.keys(): self.set_energy(kwarg["energy_eV"],"eV")
        else:
            log(logger.error,"Photon could not be initialized. It needs to be initialized with either the a given photon energy or the wavelength.")
            
    def get_energy(self,unit="J"):
        if unit == "J":
            return self._energy
        elif unit == "eV":
            return self._energy/constants.e
        else:
            log(logger.error,"%s is not a valid energy unit." % unit)

    def set_energy(self,energy,unit="J"):
        if unit == "J":
            self._energy = energy
        elif unit == "eV":
            self._energy = energy*constants.e
        else:
            log(logger.error,"%s is not a valid energy unit." % unit)

    def get_wavelength(self):
        return constants.c*constants.h/self._energy

    def set_wavelength(self,wavelength):
        self._energy = constants.c*constants.h/wavelength


        
class Profile:
    def __init__(self, model, focus_diameter):
        self.set_model(model)
        self.focus_diameter = focus_diameter
        
    def set_model(self,model):
        if model is None or model in ["top_hat","pseudo_lorentzian","gaussian"]:
            self._model = model
        else:
            log(logger.error,"Pulse profile model %s is not implemented. Change your configuration and try again.")
            exit(0)

    def get_model(self):
        return self._model

    def get_radial(self):
        if self._model is None:
            # we always hit with full power
            p = lambda r: 1. / (numpy.pi * (self.focus_diameter / 2.)**2)
        elif self._model == "top_hat":
            # focus diameter is diameter of top hat profile
            p = lambda r: (1.  / (numpy.pi * (self.focus_diameter / 2.)**2)) if r < (self.focus_diameter / 2.) else 0.
        elif self._model == "pseudo_lorentzian":
            # focus diameter is FWHM of lorentzian
            sigma = self.focus_diameter / 2.
            p = lambda r: _pseudo_lorentzian_2dnorm(r, sigma)
        elif self._model == "gaussian":
            # focus diameter is FWHM of gaussian
            sigma = self.focus_diameter / (2.*numpy.sqrt(2.*numpy.log(2.)))
            p = lambda r: _gaussian_2dnorm(r, sigma)
    
_gaussian = lambda x, sigma: numpy.exp(-x**2/(2*sigma**2))

_gaussian_2dnorm = lambda x, sigma: _gaussian(x, sigma) / ( 2 * numpy.pi * sigma**2 )

_lorentzian = lambda x, sigma: sigma**2 / (x**2 + sigma**2)

_pseudo_lorenzian_A1 = 0.74447313315648778 
_pseudo_lorenzian_A2 = 0.22788162774723308
_pseudo_lorenzian_s1 = 0.73985516665883544
_pseudo_lorenzian_s2 = 2.5588165723260907
_pseudo_lorentzian = lambda x, sigma: _pseudo_lorenzian_A1 * gaussian(x, _pseudo_lorenzian_s1*sigma) + \
                                      _pseudo_lorenzian_A2 * gaussian(x, _pseudo_lorenzian_s2*sigma)

_pseudo_lorentzian_2dnorm = lambda x, sigma: _pseudo_lorentzian(x, sigma) / ( 2. * numpy.pi * ( _pseudo_lorenzian_A1 * (_pseudo_lorenzian_s1*sigma)**2 + \
                                                                                                _pseudo_lorenzian_A2 * (_pseudo_lorenzian_s2*sigma)**2 ) )

