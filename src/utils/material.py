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
import sys, os, numpy
from scipy import constants

# Logging
import logging
logger = logging.getLogger(__name__)

# Condor
import condor
import condor._load_data
from condor.utils.log import log_and_raise_error,log_warning,log_info,log_debug

_data_dir = os.path.dirname(os.path.realpath(__file__)) + "/../data"

_atomic_scattering_factors = condor._load_data.load_atomic_scattering_factors(_data_dir)
get_atomic_scattering_factors = lambda element: _atomic_scattering_factors[element]
"""
Returns 2-dim. array of photon energy [eV] vs. real and imaginary part of the atomic scattering factor (forward scattering) for a given element.

  Args:
    :element (str): Element name (abbreviation of the latin name, for example \'He\' for helium).
"""

_atomic_masses             = condor._load_data.load_atomic_masses(_data_dir)
get_atomic_mass = lambda element: _atomic_masses[element]
"""
Returns the atomic mass (standard atomic weight in unit Dalton) for a given element.

  Args:
    :element (str): Element name (abbreviation of the latin name, for example \'He\' for helium).
"""

_atomic_numbers             = condor._load_data.load_atomic_numbers(_data_dir)
get_atomic_number = lambda element: _atomic_numbers[element]
"""
Returns the atomic number for a given element.

  Args:
    :element (str): Element name (abbreviation of the latin name, for example \'He\' for helium).
"""

atomic_names = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na''Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cp', 'Uut', 'Uuq', 'Uup', 'Uuh', 'Uus', 'Uuo']
"""
List of atom names (i.e. latin abbreviations) for all elements sorted by atomic number (increasing order).
"""


atomic_numbers = range(1,len(atomic_names)+1)
"""
List of atomic numbers of all elements in increasing order.
"""



class MaterialType:
    r"""
    Standard material types:

    ================= ====================== =================================================================== ======================
    ``material_type`` :math:`\rho_m` [kg/m3] Atomic composition                                                  Reference
    ================= ====================== =================================================================== ======================      
    ``custom``        ``massdensity``        ``atomic_composition``                                              -
    ``'water'``       995 (25 deg. C)        :math:`H_2O`                                                        [ONeil1868]_ p. 1868
    ``'protein'``     1350                   :math:`H_{86}C_{52}N_{13}O_{15}S`                                   [Bergh2008]_
    ``'dna'``         1700                   :math:`H_{11}C_{10}N_4O_6P`                                         [Bergh2008]_
    ``'lipid'``       1000                   :math:`H_{69}C_{36}O_6P`                                            [Bergh2008]_
    ``'cell'``        1000                   :math:`H_{23}C_3NO_{10}S`                                           [Bergh2008]_
    ``'poliovirus'``  1340                   :math:`C_{332652}H_{492388}N_{98245}O_{131196}P_{7501}S_{2340}`     [Molla1991]_
    ``'styrene'``     902 (25 deg. C)        :math:`C_8H_8`                                                      [Haynes2013]_ p. 3-488
    ``'sucrose'``     1581 (17 deg. C)       :math:`C_{12}H_{22O1}`                                              [Lide1998]_ p. 3-172
    ================= ====================== =================================================================== ======================
    """
    atomic_compositions = {
        'water':       { "H" :     2., "C" :     0., "N" :     0., "O" :     1., "P" :     0., "S" :     0. }, # Water H2O
        'protein':     { "H" :    86., "C" :    52., "N" :    13., "O" :    15., "P" :     0., "S" :     1. }, # Bergh et al. 2008: H86 C52 N13 O15 S
        'dna':         { "H" :    11., "C" :    10., "N" :     4., "O" :     6., "P" :     1., "S" :     0. }, # Bergh et al. 2008: H11 C10 N4 O6 P
        'lipid':       { "H" :    69., "C" :    36., "N" :     0., "O" :     6., "P" :     1., "S" :     0. }, # Bergh et al. 2008: H69 C36 O6 P
        'cell':        { "H" :    23., "C" :     3., "N" :     1., "O" :    10., "P" :     0., "S" :     1. }, # Bergh et al. 2008: H23 C3 N O10 S
        'poliovirus':  { "H" :492388., "C" :332652., "N" : 98245., "O" :131196., "P" :  7501., "S" :  2340. }, # Molla et al. 1991: C332652 H492388 N98245 0131196 P7501 S2340
        'styrene':     { "H" :     8., "C" :     8., "N" :     0., "O" :     0., "P" :     0., "S" :     0. }, # Styrene C8H8
        'sucrose':     { "H" :    22., "C" :    12., "N" :     0., "O" :    11., "P" :     0., "S" :     0. }, # Sucrose C12H22O11
    }
    """
    Dictionary of atomic compositions (available keys are the tabulated ``material_types``)
    """

    mass_densities      =  {
        'water':      995., # at 25 C: O'Neil, M.J. (ed.). The Merck Index - An Encyclopedia of Chemicals, Drugs, and Biologicals. Cambridge, UK: Royal Society of Chemistry, 2013., p. 1868
        'protein':   1350., # Bergh et al. 2008
        'dna':       1700., # Bergh et al. 2008
        'lipid':     1000., # Bergh et al. 2008
        'cell':      1000., # Bergh et al. 2008
        'poliovirus':1340., # Dans et al. 1966
        'styrene':    902., # at 25 C: Haynes, W.M. (ed.). CRC Handbook of Chemistry and Physics. 94th Edition. CRC Press LLC, Boca Raton: FL 2013-2014, p. 3-488
        'sucrose':   1581., # at 17 C: Lide, D.R. (ed.). CRC Handbook of Chemistry and Physics. 79th ed. Boca Raton, FL: CRC Press Inc., 1998-1999., p. 3-172
    }
    """
    Dictionary of mass densities (available keys are the tabulated ``material_types``)
    """

class AbstractMaterial:
    def __init__(self):
        pass

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
        scatterer_density = self.get_scatterer_density()
        
        r_0 = constants.value("classical electron radius")

        n = 1 - r_0/2/numpy.pi * photon_wavelength**2 * f * scatterer_density

        return n

    def get_transmission(self,thickness,photon_wavelength):
        r"""
        Return transmission coefficient :math:`T` for given material thickness :math:`t` and wavelength :math:`\lambda` [Henke1993]_

        .. math::

          T = e^{-\rho\,\mu_a(\lambda)\,t}

        :math:`\rho`: Average atom density

        :math:`\mu_a(\lambda)`: Photoabsorption cross section at photon energy :math:`\lambda`
        
        Args:
          :thickness (float): Material thickness in unit meter
        
          :photon_wavelength (float): Photon wavelength in unit meter

        .. [Henke1993] B.L. Henke, E.M. Gullikson, and J.C. Davis. X-ray interactions: photoabsorption, scattering, transmission, and reflection at E=50-30000 eV, Z=1-92, Atomic Data and Nuclear Data Tables Vol. 54 (no.2), 181-342 (July 1993).
          
          See also `http://henke.lbl.gov/ <http://henke.lbl.gov/>`_.
        """
        mu = self.get_photoabsorption_cross_section(photon_wavelength=photon_wavelength)
        rho = self.get_scatterer_density()

        return numpy.exp(-rho*mu*thickness)

    def get_attenuation_length(self, photon_wavelength):
        r"""
        Return the absorption length in unit meter for the given wavelength :math:`\lambda`

        .. math::
        
          \mu = \frac{1}{\rho\,\mu_a(\lambda)}
        
        :math:`\rho`: Average atom density

        :math:`\mu_a(\lambda)`: Photoabsorption cross section at photon energy :math:`\lambda`
        
        Args:
          :photon_wavelength (float): Photon wavelength in unit meter

        .. [Henke1993] B.L. Henke, E.M. Gullikson, and J.C. Davis. X-ray interactions: photoabsorption, scattering, transmission, and reflection at E=50-30000 eV, Z=1-92, Atomic Data and Nuclear Data Tables Vol. 54 (no.2), 181-342 (July 1993).
          
          See also `http://henke.lbl.gov/ <http://henke.lbl.gov/>`_.
        """
        mu = self.get_photoabsorption_cross_section(photon_wavelength=photon_wavelength)
        rho = self.get_scatterer_density()

        return (1./(mu*rho))


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

    
class ElectronDensityMaterial(AbstractMaterial):
    r"""
    Class for electron density material model
    
    Thomson scattering with the given value for the electron density is used to determine the material's scattering properties.

    Args:
      :electron_density: (float): Electron density in unit inverse cubic meter
    """
    def __init__(self, electron_density):
        AbstractMaterial.__init__(self)
        self.electron_density = electron_density

    def get_conf(self):
        conf = {}
        conf["electron_density"] = self.electron_density
        return conf
        
    def get_f(self, photon_energy):
        return complex(1., 0.)

    def get_scatterer_density(self):
        return self.electron_density
    
    
class AtomDensityMaterial(AbstractMaterial):
    r"""
    Class for material model
    
    Args:
      :material_type (str): The material type can be either ``custom`` or one of the standard types, i.e. tabulated combinations of massdensity and atomic composition, listed here :class:`condor.utils.material.MaterialType`.

    Kwargs:
      :massdensity (float): Mass density in unit kilogram per cubic meter (default ``None``)
    
      :atomic_composition (dict): Dictionary of key-value pairs for atom species (e.g. ``'H'`` for hydrogen) and concentration (default ``None``)    
    """
    def __init__(self, material_type, massdensity = None, atomic_composition = None):
        AbstractMaterial.__init__(self)
        
        self.clear_atomic_composition()

        if atomic_composition is not None and massdensity is not None and (material_type is None or material_type == "custom"):
            for element,concentration in atomic_composition.items():
                self.set_atomic_concentration(element, concentration)
            self.massdensity = massdensity

        elif material_type is not None and atomic_composition is None and massdensity is None:
            for element, concentration in MaterialType.atomic_compositions[material_type].items():
                self.set_atomic_concentration(element, concentration)
            self.massdensity = MaterialType.mass_densities[material_type]

        else:
            log_and_raise_error(logger, "Invalid arguments in Material initialization.")

    def get_conf(self):
        conf = {}
        conf["material_type"]      = "custom"
        conf["atomic_composition"] = self.get_atomic_composition()
        conf["massdensity"]        = self.massdensity
        return conf

    def clear_atomic_composition(self):
        """
        Empty atomic composition dictionary
        """
        self._atomic_composition = {}
    
    def set_atomic_concentration(self, element, relative_concentration):
        r"""
        Set the concentration of a given atomic species

        Args:
          :element (str): Atomic species (e.g. ``'H'`` for hydrogen)

          :relative_concentration (float): Relative quantity of atoms of the given atomic species with respect to the others (e.g. for water: hydrogen concentration ``2.``, oxygen concentration ``1.``)
        """
        if element not in atomic_names:
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
        
    def get_f(self, photon_wavelength):
        r"""
        Get effective average complex scattering factor for forward scattering at a given photon wavlength from Henke tables

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

    def get_scatterer_density(self):
        r"""
        Return total atom density :math:`\rho` in unit inverse cubic meter

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
            # sum up average mass
            M += atomic_composition[element]*get_atomic_mass(element)*u

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
            M += atomic_composition[element]*get_atomic_mass(element)*u
            Q += atomic_composition[element]*get_atomic_number(element)

        electron_density = Q*self.massdensity/M
        
        return electron_density


def get_f_element(element, photon_energy_eV):
    r"""
    Get the scattering factor for an element through linear interpolation of the tabulated values (Henke tables)

    Args:
      :element (str): Atomic species (e.g. ``'H'`` for hydrogen)

      :photon_energy_eV: Photon energy in unit eV
    """
    
    SF_X = get_atomic_scattering_factors(element)
    f1 = numpy.interp(photon_energy_eV,SF_X[:,0],SF_X[:,1])
    f2 = numpy.interp(photon_energy_eV,SF_X[:,0],SF_X[:,2])

    return complex(f1,f2)


class MaterialMap:
    def __init__(self, shape):
        if len(shape) != 3:
            log_and_raise_error(logger, "%s is an invald shape for initialisation of MaterialMap.", str(shape))
        self._shape = tuple(shape)

    def add_material(self, material, density_map):
        if not isinstance(material, Material):
            log_and_raise_error(logger, "Cannot add material %s. It is not an instance of Material." % str(material))
        if density_map.shape != self._shape:
            log_and_raise_error(logger, "Cannot add material. Density map has incompatible shape: %s. Should be %s." % (str(density_map.shape), str(self._shape)))
        self.materials.append(material)
        self.density_maps.append(density_map)

    def get_n(self, photon_wavelength):
        dn = self.get_dn(photon_wavelength)
        n = 1 - dn
        return n
        
    def get_dn(self, photon_wavelength):
        dn = numpy.zeros(shape=self._shape, dtype=numpy.complex128)
        for mat,dmap in zip(self.materials, self.density_maps):
            dn += mat.get_dn(photon_wavelength) * dmap
        return dn

    def get_beta(self, photon_wavelength):
        dn = self.get_dn(photon_wavelength)
        return dn.imag

    def get_delta(self, photon_wavelength):
        dn = self.get_dn(photon_wavelength)
        return dn.real
    
    def get_photoabsorption_cross_section(self, photon_wavelength):
        pacs = numpy.zeros(shape=self._shape, dtype=numpy.float64)
        for mat,dmap in zip(self.materials, self.density_maps):
            pacs += mat.get_photoabsorption_cross_section(photon_wavelength) * dmap
        return pacs
                
    def get_f(self, photon_wavelength):
        f = numpy.zeros(shape=self._shape, dtype=numpy.complex128)
        for mat,dmap in zip(self.materials, self.density_maps):
            f += mat.get_f(photon_wavelength)*dmap
        return trans
        
    def get_electron_density(self, photon_wavelength):
        ed = numpy.zeros(shape=self._shape, dtype=numpy.float64)
        for mat,dmap in zip(self.materials, self.density_maps):
            ed += mat.get_electron_density(photon_wavelength) * dmap
        return ed


#class DensityMap:
#    
#    def __init__(self, shape):
#        self.density = numpy.zeros(shape=(shape[0], shape[1], shape[2], len(atomic_numbers.keys())),dtype=numpy.float64)
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
#        for (el, de) in zip(atomic_numbers.keys(), self.density):
#            s += de * get_f_element(el, photon_energy_eV)
#
#        n = 1 - r_0 / (2*numpy.pi) * wavelength**2 * s
#
#        return n
#
#    def get_dn(self, wavelength):
#        return (1-self.get_n(wavelength))
