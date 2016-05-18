import unittest
from condor.utils import material
import numpy
import scipy.constants as constants

class TestCaseMaterial(unittest.TestCase):
    def test_henke_constants(self):
        # Output from http://henke.lbl.gov/ for He at 1500 eV
        energy_eV = 1500.
        element = 'He'
        # Atomic Properties
        # Atomic Weight : 4.003
        # Photoabsorption Cross Section : 16.31 cm^2/g
        # Inelastic Cross Section : 0.2135E-01 cm^2/g
        # f1 : 2.006
        f1_expected = 2.006
        # f2 : 0.2327E-02
        f2_expected = 0.2327E-02
        # Material Properties
        # Density : 0.179E-03 g/cm^3
        # Delta : 1.6511E-08
        delta_expected = 1.6511E-08
        # Beta : 1.9146E-11
        beta_expected = 1.9146E-11
        # Attenuation Length : 3.4311E+06 microns
        mu_expected = 3.4311E+06 * 1E-6
        arr_energy_eV = material.get_atomic_scattering_factors(element)[:,0]
        arr_f1 = material.get_atomic_scattering_factors('He')[:,1]
        arr_f2 = material.get_atomic_scattering_factors('He')[:,2]
        f1 = numpy.interp(energy_eV, arr_energy_eV, arr_f1)
        f2 = numpy.interp(energy_eV, arr_energy_eV, arr_f2)
        self.assertAlmostEqual(f1 , f1_expected, 3)
        self.assertAlmostEqual(f2/1E-2, f2_expected/1E-2, 4)
        f1 = material.get_f_element(element, energy_eV).real
        f2 = material.get_f_element(element, energy_eV).imag
        self.assertAlmostEqual(f1 , f1_expected, 3)
        self.assertAlmostEqual(f2/1E-2, f2_expected/1E-2, 4)
        

    def test_material(self):
        # Output from http://henke.lbl.gov/ for He at 1500 eV
        energy_eV = 1500.
        wavelength = constants.c*constants.h/(energy_eV*constants.e)
        element = 'He'
        # Atomic Properties
        # Atomic Weight : 4.003
        # Photoabsorption Cross Section : 16.31 cm^2/g
        # Inelastic Cross Section : 0.2135E-01 cm^2/g
        # f1 : 2.006
        f1_expected = 2.006
        # f2 : 0.2327E-02
        f2_expected = 0.2327E-02
        # Material Properties
        # Density : 0.179E-03 g/cm^3
        massdensity = 0.179E-03 * 1E-3/(1E-2)**3
        # Delta : 1.6511E-08
        delta_expected = 1.6511E-08
        # Beta : 1.9146E-11
        beta_expected = 1.9146E-11
        # Attenuation Length : 3.4311E+06 microns
        mu_expected = 3.4311E+06 * 1E-6
        M = material.AtomDensityMaterial(material_type='custom', massdensity=massdensity, atomic_composition={'He' : 1})
        f1 = M.get_f(wavelength).real
        f2 = M.get_f(wavelength).imag
        delta = M.get_delta(wavelength)
        beta = M.get_beta(wavelength)
        mu = M.get_attenuation_length(wavelength)
        self.assertAlmostEqual(f1 , f1_expected, 3)
        self.assertAlmostEqual(f2/1E-2 , f2_expected/1E-2, 4)
        self.assertAlmostEqual(delta/1E-8 , delta_expected/1E-8, 2)
        self.assertAlmostEqual(mu , mu_expected, 1)
        
