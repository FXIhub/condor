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

import sys,numpy
import logging
logger = logging.getLogger("Condor")
if "utils" not in sys.path: sys.path.append("utils")
from scipy import constants
import config

class Material:
    """
    A class of a sample object.
    Sample material.

    """
    def __init__(self,**kwargs):
        # Check for valid keyword arguments
        allk = ["massdensity","material_type"]
        self._unproc_kws = [k for k in kwargs.keys() if (k not in allk and k[0] != "c" and len(k) > 3)]
        if len(self._unproc_kws) > 0:
            print self._unproc_kws
            logger.error("Material object initialisation failed due to illegal keyword arguments.")
            return
        # Start intialisation
        if "massdensity" in kwargs:
            self.materialtype = "custom"
            for key in kwargs:
                if key[0] == 'c' or key == 'massdensity':
                    exec "self." + key + " = kwargs[key]"
                else:
                    logger.error("%s is no valid argument for custom initialization of Material." % key)
                    return
        elif "material_type" in kwargs:
            self.material_type = kwargs['material_type']
            self.massdensity = config.DICT_massdensity[self.material_type]
            self.cH = config.DICT_atomic_composition[self.material_type][0]
            self.cC = config.DICT_atomic_composition[self.material_type][1]
            self.cN = config.DICT_atomic_composition[self.material_type][2]
            self.cO = config.DICT_atomic_composition[self.material_type][3]
            self.cP = config.DICT_atomic_composition[self.material_type][4]
            self.cS = config.DICT_atomic_composition[self.material_type][5]
            self.cAu = config.DICT_atomic_composition[self.material_type][6]
        else:
            logger.error("No valid arguments for Material initialization.")
            return

 
    def get_n(self,photon_wavelength):
        """
        Obtains complex refractive index.
        Henke (1994): n = 1 - r_0/(2pi) lambda^2 sum_q rho_q f_q(0)
        r_0: classical electron radius
        rho_q: atomic number density of atom species q
        f_q(0): atomic scattering factor (forward scattering) of atom species q
        """

        f = self.get_f(photon_wavelength)
        atom_density = self.get_atom_density()
        
        r_0 = constants.value("classical electron radius")

        n = 1 - r_0/2/numpy.pi * photon_wavelength**2 * f * atom_density

        return n

    def get_dn(self,photon_wavelength=None):
        return (1-self.get_n(photon_wavelength))

    # convenience functions
    # n = 1 - delta - i beta
    def get_delta(self,photon_wavelength):
        return (1-self.get_n(photon_wavelength=photon_wavelength).real)
    def get_beta(self,photon_wavelength):
        return (-self.get_n(photon_wavelength=photon_wavelength).imag)

    def get_photoabsorption_cross_section(self,photon_wavelength):

        r_0 = constants.value("classical electron radius")
        h =  constants.h
        c =  constants.c
        qe = constants.e

        mu = 2*r_0*photon_wavelength*self.get_f(photon_wavelength).imag

        return mu

    def get_transmission(self,thickness,photon_wavelength):

        n = self.get_n(photon_wavelength)
        mu = self.get_photoabsorption_cross_section(photon_wavelength=photon_wavelength)
        rho = self.get_atom_density()

        return numpy.exp(-rho*mu*thickness)

    def get_f(self,photon_wavelength):
    
        atomic_composition = self.get_atomic_composition_dict()

        r_0 = constants.value("classical electron radius")
        h   =  constants.h
        c   =  constants.c
        qe   = constants.e
        photon_energy_eV = h*c/photon_wavelength/qe
        
        f_sum = complex(0.,0.)
        for element in atomic_composition.keys():
            # sum up average atom factor
            f = get_f_element(element,photon_energy_eV)
            f_sum += atomic_composition[element] * f
        
        return f_sum


    def get_atom_density(self):
                
        u = constants.value("atomic mass constant")

        atomic_composition = self.get_atomic_composition_dict()

        M = 0
        for element in atomic_composition.keys():
            # sum up mass
            M += atomic_composition[element]*config.DICT_atomic_mass[element]*u

        number_density = self.massdensity/M
        
        return number_density


    def get_electron_density(self):

        u = constants.value("atomic mass constant")

        atomic_composition = self.get_atomic_composition_dict()

        M = 0
        Q = 0
        for element in atomic_composition.keys():
            # sum up electrons
            M += atomic_composition[element]*config.DICT_atomic_mass[element]*u
            Q += atomic_composition[element]*config.DICT_atomic_number[element]

        electron_density = Q*self.massdensity/M
        
        return electron_density
        
        
    def get_atomic_composition_dict(self):

        atomic_composition = {}
        
        for key in self.__dict__.keys():
            if key[0] == 'c' and key[1:] in config.DICT_atomic_number.keys():
                exec "c_tmp = self." + key
                atomic_composition[key[1:]] = c_tmp 
 
        tmp_sum = float(sum(atomic_composition.values()))
        for element in atomic_composition.keys():
            atomic_composition[element] /= tmp_sum 
        
        return atomic_composition



    

def get_f_element(element, photon_energy_eV):
    """
    Get the scattering factor for an element through linear interpolation.
    """
    
    SF_X = config.DICT_scattering_factors[element]
    f1 = numpy.interp(photon_energy_eV,SF_X[:,0],SF_X[:,1])
    f2 = numpy.interp(photon_energy_eV,SF_X[:,0],SF_X[:,2])

    return complex(f1,f2)


class DensityMap:
    
    def __init__(self, shape):
        self.density = numpy.zeros(shape=(shape[0], shape[1], shape[2], len(config.DICT_atomic_number.keys())),dtype=numpy.float64)

    def get_n(self, wavelength):
        """
        Obtains complex refractive index.
        Henke (1994): n = 1 - r_0/(2pi) lambda^2 sum_q rho_q f_q(0)
        r_0: classical electron radius
        rho_q: atomic number density of atom species q
        f_q(0): atomic scattering factor (forward scattering) of atom species q
        """

        r_0 = constants.value("classical electron radius")
        h   =  constants.h
        c   =  constants.c
        qe   = constants.e
        photon_energy_eV = h*c/photon_wavelength/qe

        s = numpy.zeros(shape=(shape[0], shape[1], shape[2]), dtype=numpy.complex128)
        for (el, de) in zip(config.DICT_atomic_number.keys(), self.density):
            s += de * get_f_element(el, photon_energy_eV)

        n = 1 - r_0 / (2*numpy.pi) * wavelength**2 * s

        return n

    def get_dn(self, wavelength):
        return (1-self.get_n(wavelength))
