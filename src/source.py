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

import numpy
import config,logging
logger = logging.getLogger('Condor')
from scipy import constants

from variation import Variation

class Source:
    """
    A subclass of the input object.
    Defines properties of the FEL x-ray pulse.

    """
    def __init__(self,**kwargs):       
        # Check for required keyword arguments
        reqk = ["wavelength","focus_diameter","pulse_energy"]
        for k in reqk:
            if k not in kwargs.keys():
                logger.error("Cannot initialize Source instance. %s is a necessary keyword." % k)
                return
        # Check for valid keyword arguments
        allk = ["parent",
                "wavelength","focus_diameter","pulse_energy","pulse_energy_mean",
                "pulse_energy_variation","pulse_energy_spread","pulse_energy_variation_n"]
        self._unproc_kws = [k for k in kwargs.keys() if k not in allk]
        if len(self._unproc_kws) > 0:
            print self._unproc_kws
            logger.error("Source object initialisation failed due to illegal keyword arguments.")
            return
        # Start initialisation
        self._parent = kwargs.get('parent',None)
        self.photon = Photon(wavelength=kwargs["wavelength"])
        self.focus_diameter = kwargs["focus_diameter"]
        if "pulse_energy" in kwargs:
            # Maintain depreciated keyword
            self.pulse_energy_mean = kwargs["pulse_energy"]
        else:
            self.pulse_energy_mean = kwargs["pulse_energy_mean"]
        self.set_pulse_energy_variation(variation=kwargs.get("pulse_energy_variation",None),
                                        spread=kwargs.get("pulse_energy_spread",None),
                                        n=kwargs.get("pulse_energy_variation_n",None))
        logger.debug("Source configured")

    def set_pulse_energy_variation(self,variation=None, spread=None, n=None):
        self._pulse_energy_variation = Variation(variation,spread,n,number_of_dimensions=1,name="pulse energy")

    def get_intensity(self,unit="ph/m2"):
        I = self.pulse_energy / self.get_area()
        if unit == "J/m2":
            pass
        elif unit == "ph/m2":
            I /= self.photon.get_energy("J") 
        elif unit == "J/um2":
            I *= 1.E-12
        elif unit == "mJ/um2":
            I *= 1.E-9
        else:
            logger.error("%s is not a valid unit." % unit)
            return
        return I

    def get_area(self):
        return numpy.pi*(self.focus_diameter/2.0)**2

    def get_next(self):
        self._next_pulse_energy()
        return {"pulse_energy":self.pulse_energy,
                "wavelength":self.photon.get_wavelength(),
                "intensity":self.get_intensity(),
                "intensity_ph_m2":self.get_intensity("ph/m2")}

    def _next_pulse_energy(self):
        p = self._pulse_energy_variation.get(self.pulse_energy_mean)
        # Non-random
        if self._pulse_energy_variation._mode in [None,"range"]:
            if p <= 0:
                logger.error("Pulse energy smaller-equals zero. Change your configuration.")
            else:
                self.pulse_energy = p
        # Random
        else:
            if p <= 0.:
                logger.warning("Pulse energy smaller-equals zero. Try again.")
                self._next_pulse_energy()
            else:
                self.pulse_energy = p


    
class Photon:
    def __init__(self,**kwarg):
        if "wavelength" in kwarg.keys(): self.set_wavelength(kwarg["wavelength"])
        elif "energy" in kwarg.keys(): self.set_energy(kwarg["energy"],"J")
        elif "energy_eV" in kwarg.keys(): self.set_energy(kwarg["energy_eV"],"eV")
        else:
            logger.error("Photon could not be initialized. It needs to be initialized with either the a given photon energy or the wavelength.")
            
    def get_energy(self,unit="J"):
        if unit == "J":
            return self._energy
        elif unit == "eV":
            return self._energy/constants.e
        else:
            logger.error("%s is not a valid energy unit." % unit)

    def set_energy(self,energy,unit="J"):
        if unit == "J":
            self._energy = energy
        elif unit == "eV":
            self._energy = energy*constants.e
        else:
            logger.error("%s is not a valid energy unit." % unit)

    def get_wavelength(self):
        return constants.c*constants.h/self._energy

    def set_wavelength(self,wavelength):
        self._energy = constants.c*constants.h/wavelength


            


            
