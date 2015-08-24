# -----------------------------------------------------------------------------------------------------
# CONOR
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

# System packages
import sys, numpy

# Logging
import logging
logger = logging.getLogger(__name__)
import condor
import condor.utils.log
from condor.utils.log import log_and_raise_error,log_warning,log_info,log_debug

# Constants
from scipy import constants

class Material:
    def __init__(self, material_type = None, massdensity = None, atomic_composition = None):

        self.clear_atomic_composition()
        
        if atomic_composition is not None and massdensity is not None and (material_type is None or material_type == "custom"):
            for element,quantity in atomic_composition.items():
                self.add_atomic_species(element, quantity)
            self.massdensity = massdensity

        elif material_type is not None and atomic_composition is None and massdensity is None:
            for element,quantity in condor.CONDOR_atomic_compositions[material_type].items():
                self.add_atomic_species(element,quantity)
            self.massdensity = condor.CONDOR_mass_densities[material_type]

        else:
            log_and_raise_error(logger, "Invalid arguments in Material initialization.")

    def clear_atomic_composition(self):
        self._atomic_composition = {}
            
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
    
        atomic_composition = self.get_atomic_composition(normed=True)

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

        atomic_composition = self.get_atomic_composition(normed=True)

        M = 0
        for element in atomic_composition.keys():
            # sum up mass
            M += atomic_composition[element]*condor.CONDOR_atomic_masses[element]*u

        number_density = self.massdensity/M
        
        return number_density


    def get_electron_density(self):

        u = constants.value("atomic mass constant")

        atomic_composition = self.get_atomic_composition(normed=True)

        M = 0
        Q = 0
        for element in atomic_composition.keys():
            # sum up electrons
            M += atomic_composition[element]*condor.CONDOR_atomic_masses[element]*u
            Q += atomic_composition[element]*condor.CONDOR_atomic_numbers[element]

        electron_density = Q*self.massdensity/M
        
        return electron_density
        

    def add_atomic_species(self, element, quantity):
        if element not in condor.CONDOR_atomic_numbers:
            log_and_raise_error(logger, "Cannot add element \"%s\". Invalid name." % element)
        self._atomic_composition[element] = quantity
    
    def get_atomic_composition(self, normed=False):
        
        atomic_composition = self._atomic_composition.copy()

        if normed:
            s = numpy.array(atomic_composition.values(), dtype=numpy.float64).sum()
            for element in atomic_composition.keys():
                atomic_composition[element] /= s 

        return atomic_composition



    

def get_f_element(element, photon_energy_eV):
    """
    Get the scattering factor for an element through linear interpolation.
    """
    
    SF_X = condor.CONDOR_atomic_scattering_factors[element]
    f1 = numpy.interp(photon_energy_eV,SF_X[:,0],SF_X[:,1])
    f2 = numpy.interp(photon_energy_eV,SF_X[:,0],SF_X[:,2])

    return complex(f1,f2)


class DensityMap:
    
    def __init__(self, shape):
        self.density = numpy.zeros(shape=(shape[0], shape[1], shape[2], len(condor.CONDOR_atomic_numbers.keys())),dtype=numpy.float64)

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
        for (el, de) in zip(condor.CONDOR_atomic_numbers.keys(), self.density):
            s += de * get_f_element(el, photon_energy_eV)

        n = 1 - r_0 / (2*numpy.pi) * wavelength**2 * s

        return n

    def get_dn(self, wavelength):
        return (1-self.get_n(wavelength))
