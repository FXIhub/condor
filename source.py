import pylab
import config

class Source:
    """
    A subclass of the input object.
    Defines properties of the FEL x-ray pulse.

    """
    def __init__(self,**kwargs):
        self._parent = kwargs.get('parent',None)
        self.photon = Photon(wavelength=kwargs.get('wavelength',5.7E-09))
        self.focus_diameter = kwargs.get('focus_diameter',20E-06)
        self.pulse_energy = kwargs.get('pulse_energy',100E-06)

    def get_area(self):
        return pylab.pi*(self.focus_diameter/2.0)**2
    
class Photon:
    def __init__(self,**kwarg):
        if "wavelength" in kwarg.keys(): self.set_wavelength(kwarg["wavelength"])
        elif "energy" in kwarg.keys(): self.set_energy(kwarg["energy"])
        elif "energy_eV" in kwarg.keys(): self.set_energy_eV(kwarg["energy_eV"])
        else:
            print "ERROR: Photon needs to be initialized with either the energy or the wavelength."
            
    def get_energy(self,unit="J"):
        if unit == "J":
            return self._energy
        elif unit == "eV":
            return self._energy/config.DICT_physical_constants["e"]
        else:
            print "ERROR: %s is not a valid energy unit." % unit

    def set_energy(self,energy,unit="J"):
        if unit == "J":
            self._energy = energy
        elif unit == "eV":
            self._energy = energy*config.DICT_physical_constants["e"]
        else:
            print "ERROR: %s is not a valid energy unit." % unit

    def get_wavelength(self):
        return config.DICT_physical_constants["c"]*config.DICT_physical_constants["h"]/self._energy

    def set_wavelength(self,wavelength):
        self._energy = config.DICT_physical_constants["c"]*config.DICT_physical_constants["h"]/wavelength


            


            
