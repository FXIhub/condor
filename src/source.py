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

class Source:
    """
    A subclass of the input object.
    Defines properties of the FEL x-ray pulse.

    """
    def __init__(self,**kwargs):
        self._parent = kwargs.get('parent',None)
        reqk = ["wavelength","focus_diameter","pulse_energy"]
        for k in reqk:
            if k not in kwargs.keys():
                logger.error("Cannot initialize Source instance. %s is a necessary keyword." % k)
                return
        self.photon = Photon(wavelength=kwargs["wavelength"])
        self.focus_diameter = kwargs["focus_diameter"]
        if "pulse_energy" in kwargs:
            # Maintain depreciated keyword
            self.pulse_energy_mean = kwargs["pulse_energy"]
        else:
            self.pulse_energy_mean = kwargs["pulse_energy_mean"]
        self.pulse_energy_variation = kwargs.get("pulse_energy_variation",None)
        self.pulse_energy_spread = kwargs.get("pulse_energy_spread",0.)
        self._next()
        logger.debug("Source configured")

    def get_intensity(self,unit="ph/m2"):
        I = self.pulse_energy / self.get_area()
        if unit == "J/m2":
            pass
        if unit == "ph/m2":
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

    def _next(self):
        if self.pulse_energy_variation is None:
            self.pulse_energy = self.pulse_energy_mean
        else:
            if self.pulse_energy_variation == "normal":
                self.pulse_energy = numpy.random.normal(self.pulse_energy_mean,self.pulse_energy_spread)
            elif self.pulse_energy_variation == "uniform":
                self.pulse_energy = numpy.random.uniform(self.pulse_energy_mean-self.pulse_energy_spread/2.,self.pulse_energy_mean+self.pulse_energy_spread/2.)
            else:
                logger.error("%s is an illegal variation for the pulse energy.",self.pulse_energy_variation)
            if self.pulse_energy <= 0.:
                self._next()

    
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


            


            
