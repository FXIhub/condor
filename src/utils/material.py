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
    r"""
    Class for material model
    
    **Arguments:**

      :material_type (str): The material type can be either ``custom`` or one of the standard types, i.e. tabulated combinations of massdensity and atomic composition, listed here:

        ================= ====================== =================================================================== ===================
        ``material_type`` :math:`\rho_m` [kg/m3] Atomic composition                                                  Reference
        ================= ====================== =================================================================== ===================        
        ``custom``        ``massdensity``        ``atomic_composition``                                              -
        ``'water'``       995 (25 deg. C)        :math:`H_2O`                                                        [ONeil1868]_
        ``'protein'``     1350                   :math:`H_{86}C_{52}N_{13}O_{15}S`                                   [Bergh2008]_
        ``'dna'``         1700                   :math:`H_{11}C_{10}N_4O_6P`                                         [Bergh2008]_
        ``'lipid'``       1000                   :math:`H_{69}C_{36}O_6P`                                            [Bergh2008]_
        ``'cell'``        1000                   :math:`H_{23}C_3NO_{10}S`                                           [Bergh2008]_
        ``'poliovirus'``  1340                   :math:`C_{332652}H_{492388}N_{98245}O_{131196}P_{7501}S_{2340}`     [Molla1991]_
        ``'styrene'``     902 (25 deg. C)        :math:`C_8H_8`                                                      [Haynes2013]_
        ``'sucrose'``     1581 (17 deg. C)       :math:`C_{12}H_{22O1}`                                              [Lide1998]_
        ================= ====================== =================================================================== ===================

    **Keyword arguments:**

      :massdensity (float): Mass density in unit kilogram per cubic meter (default ``None``)
    
      :atomic_composition (dict): Dictionary of key-value pairs for atom species (e.g. ``'H'`` for hydrogen) and concentration (default ``None``)    

    .. [ONeil1868] O'Neil, M.J. (ed.). The Merck Index - An Encyclopedia of Chemicals, Drugs, and Biologicals. Cambridge, UK: Royal Society of Chemistry, 2013., p. 1868
    .. [Bergh2008] Bergh et al. 2008
    .. [Molla1991] Molla et al. 1991
    .. [Haynes2013] Haynes, W.M. (ed.). CRC Handbook of Chemistry and Physics. 94th Edition. CRC Press LLC, Boca Raton: FL 2013-2014, p. 3-488
    .. [Lide1998] Lide, D.R. (ed.). CRC Handbook of Chemistry and Physics. 79th ed. Boca Raton, FL: CRC Press Inc., 1998-1999., p. 3-172
    """
    def __init__(self, material_type, massdensity = None, atomic_composition = None):

        self.clear_atomic_composition()
        
        if atomic_composition is not None and massdensity is not None and (material_type is None or material_type == "custom"):
            for element,concentration in atomic_composition.items():
                self.set_atomic_concentration(element, concentration)
            self.massdensity = massdensity

        elif material_type is not None and atomic_composition is None and massdensity is None:
            for element,quantity in condor.CONDOR_atomic_compositions[material_type].items():
                self.add_atomic_species(element,quantity)
            self.massdensity = condor.CONDOR_mass_densities[material_type]

        else:
            log_and_raise_error(logger, "Invalid arguments in Material initialization.")

    def clear_atomic_composition(self):
        """
        Empty atomic composition dictionary
        """
        self._atomic_composition = {}
            
    def get_n(self,photon_wavelength):
        r"""
        Return complex refractive index at a given wavelength (Henke, 1994)

        .. math::

          n = 1 - \frac{ r_0 }{ 2\pi } \lambda^2 \sum_i \rho_i f_i(0)
        
        :math:`r_0`: classical electron radius

        :math:`\rho_q`: atomic number density of atom species :math:`i`

        :math:`f_q(0)`: atomic scattering factor (forward scattering) of atom species :math:`i`

        Args:
         
          :photon_wavelength (float): Photon wavelength in unit meter
        """
        f = self.get_f(photon_wavelength)
        atom_density = self.get_atom_density()
        
        r_0 = constants.value("classical electron radius")

        n = 1 - r_0/2/numpy.pi * photon_wavelength**2 * f * atom_density

        return n

    def get_dn(self, photon_wavelength):
        r"""
        Return :math:`\delta n` at a given wavelength

        .. math::

          \delta n = 1 - n

        :math:`n`: Refractive index
        
        See also :meth:`condor.utils.material.Material.get_n`

        Args:

          :photon_wavelength (float): Photon wavelength in unit meter

        """
        return (1-self.get_n(photon_wavelength))

    # convenience functions
    # n = 1 - delta - i beta
    def get_delta(self, photon_wavelength):
        r"""
        Return :math:`\delta` (real part of :math:`\delta n`) at a given wavelength

        .. math:

          n = 1 - \delta n = 1 - \delta - i \beta

        :math:`n`: Refractive index
        
        See also :meth:`condor.utils.material.Material.get_n`

        Args:

          :photon_wavelength (float): Photon wavelength in unit meter

        """        
        return self.get_dn(photon_wavelength=photon_wavelength).real

    def get_beta(self, photon_wavelength):
        r"""
        Return :math:`\beta` (imaginary part of :math:`\delta n`) at a given wavelength

        .. math:: 
        
          n = 1 - \delta n = 1 - \delta - i \beta

        :math:`n`: Refractive index
        
        See also :meth:`condor.utils.material.Material.get_n`

        Args:

          :photon_wavelength (float): Photon wavelength in unit meter

        """        
        return self.get_dn(photon_wavelength=photon_wavelength).imag

    def get_photoabsorption_cross_section(self, photon_wavelength):
        r"""
        Return the photoabsorption cross section :math:`\mu_a` at a given wavelength :math:`\lambda`

        .. math::
        
          \mu_a = 2 r_0 \lambda f_2

        :math:`r_0`: classical electron radius

        :math:`f_2`: imaginary part of the atomic scattering factor     

        Args:

          :photon_wavelength (float): Photon wavelength in unit meter
        """        
        r_0 = constants.value("classical electron radius")
        h =  constants.h
        c =  constants.c
        qe = constants.e

        mu = 2*r_0*photon_wavelength*self.get_f(photon_wavelength).imag

        return mu

    def get_transmission(self,thickness,photon_wavelength):
        r"""
        Return transmission coefficient :math:`T` for given material thickness :math:`t` at a given wavelength :math:`\lambda` (Henke, 1993)

        .. math::

          T = e^{-\rho\,\mu_a(\lambda)\,t}

        :math:`\rho`: Average atom density

        :math:`\mu_a(\lambda)`: Photoabsorption cross section at photon energy :math:`\lambda`
        
        Args:

          :thickness (float): Material thickness in unit meter
        
          :photon_wavelength (float): Photon wavelength in unit meter
        """

        n = self.get_n(photon_wavelength)
        mu = self.get_photoabsorption_cross_section(photon_wavelength=photon_wavelength)
        rho = self.get_atom_density()

        return numpy.exp(-rho*mu*thickness)

    def get_f(self, photon_wavelength):
        r"""
        Read complex scattering factor at a given photon wavlength from Henke tables

        Args:
    
          :photon_wavlength (float): Photon wavelength in unit meter
        """
    
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
        r"""
        Return average atom density :math:`\rho` weighted by standard atomic mass in unit inverse cubic meter

        .. math::

          \rho = \frac{\rho_m}{\sum_i c_i m_i}
        
        :math:`\rho_m`: Mass denisty of material

        :math:`c_i`: Normalised fraction of atom species :math:`i`

        :math:`m_i`: Standard atomic mass of atom species :math:`i`       
        """
                
        u = constants.value("atomic mass constant")

        atomic_composition = self.get_atomic_composition(normed=True)

        M = 0
        for element in atomic_composition.keys():
            # sum up mass
            M += atomic_composition[element]*condor.CONDOR_atomic_masses[element]*u

        number_density = self.massdensity/M
        
        return number_density


    def get_electron_density(self):
        r"""
        Return electron density :math:`\rho_e` in unit inverse cubic meter

        .. math::

          \rho_e = \frac{\rho_m \cdot \sum_i}{\left( \sum_i c_i m_i \right) \left( \sum_i c_i Z_i \right)}

        :math:`\rho_m`: Mass denisty of material

        :math:`c_i`: Normalised fraction of atom species :math:`i`

        :math:`m_i`: Standard atomic mass of atom species :math:`i` 
      
        :math:`Z_i`: Atomic number of atom species :math:`i`
        """
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
        

    def set_atomic_concentration(self, element, relative_concentration):
        r"""
        Set the concentration of a given atomic species

        Args:
        
          :element (str): Atomic species (e.g. ``'H'`` for hydrogen)

          :relative_concentration (float): Relative quantity of atoms of the given atomic species with respect to the others (e.g. for water: hydrogen concentration 2., oxygen concentration 1.)
        """
        if element not in condor.CONDOR_atomic_numbers:
            log_and_raise_error(logger, "Cannot add element \"%s\". Invalid name." % element)
        self._atomic_composition[element] = relative_concentration
    
    def get_atomic_composition(self, normed=False):
        r"""
        Return dictionary of atomic concentrations

        Args:

          :normed (bool): If ``True`` the concentrations are rescaled by a common factor such that their sum equals 1 (default ``False``)
        """
        
        atomic_composition = self._atomic_composition.copy()

        if normed:
            s = numpy.array(atomic_composition.values(), dtype=numpy.float64).sum()
            for element in atomic_composition.keys():
                atomic_composition[element] /= s 

        return atomic_composition



    

def get_f_element(element, photon_energy_eV):
    r"""
    Get the scattering factor for an element through linear interpolation of the tabulated values (Henke tables)

    Args:
    
      :element (str): Atomic species (e.g. ``'H'`` for hydrogen)

      :photon_energy_eV: Photon energy in unit eV
    """
    
    SF_X = condor.CONDOR_atomic_scattering_factors[element]
    f1 = numpy.interp(photon_energy_eV,SF_X[:,0],SF_X[:,1])
    f2 = numpy.interp(photon_energy_eV,SF_X[:,0],SF_X[:,2])

    return complex(f1,f2)


#class DensityMap:
#    
#    def __init__(self, shape):
#        self.density = numpy.zeros(shape=(shape[0], shape[1], shape[2], len(condor.CONDOR_atomic_numbers.keys())),dtype=numpy.float64)
#
#    def get_n(self, wavelength):
#        """
#        Obtains complex refractive index.
#        Henke (1994): n = 1 - r_0/(2pi) lambda^2 sum_q rho_q f_q(0)
#        r_0: classical electron radius
#        rho_q: atomic number density of atom species q
#        f_q(0): atomic scattering factor (forward scattering) of atom species q
#        """
#
#        r_0 = constants.value("classical electron radius")
#        h   =  constants.h
#        c   =  constants.c
#        qe   = constants.e
#        photon_energy_eV = h*c/photon_wavelength/qe
#
#        s = numpy.zeros(shape=(shape[0], shape[1], shape[2]), dtype=numpy.complex128)
#        for (el, de) in zip(condor.CONDOR_atomic_numbers.keys(), self.density):
#            s += de * get_f_element(el, photon_energy_eV)
#
#        n = 1 - r_0 / (2*numpy.pi) * wavelength**2 * s
#
#        return n
#
#    def get_dn(self, wavelength):
#        return (1-self.get_n(wavelength))
