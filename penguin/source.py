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
        self.pulse_energy = kwargs["pulse_energy"]
        logger.debug("Source configured")

    def get_intensity(self,unit="J/m2"):
        I = self.pulse_energy / self.get_area()
        if unit == "ph/m2":
            I /= self.photon.get_energy("J") 
        elif unit == "J/um2":
            I *= 1.E-12
        else:
            logger.error("%s is not a valid unit." % unit)
            return
        return I

    def get_area(self):
        return numpy.pi*(self.focus_diameter/2.0)**2
    
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
            return self._energy/config.DICT_physical_constants["e"]
        else:
            logger.error("%s is not a valid energy unit." % unit)

    def set_energy(self,energy,unit="J"):
        if unit == "J":
            self._energy = energy
        elif unit == "eV":
            self._energy = energy*config.DICT_physical_constants["e"]
        else:
            logger.error("%s is not a valid energy unit." % unit)

    def get_wavelength(self):
        return config.DICT_physical_constants["c"]*config.DICT_physical_constants["h"]/self._energy

    def set_wavelength(self,wavelength):
        self._energy = config.DICT_physical_constants["c"]*config.DICT_physical_constants["h"]/wavelength


            


            
