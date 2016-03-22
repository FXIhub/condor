import unittest
from condor.utils import photon

class TestCasePhoton(unittest.TestCase):
    def test_consistency(self):
        wavelength = 1E-10
        energy_eV  = 12.39841662513396E3
        energy     = energy_eV*1.60217662E-19
        frequency  = energy/6.62607004E-34
        for P in [photon.Photon(wavelength=wavelength), photon.Photon(energy=energy), photon.Photon(energy_eV=energy_eV), photon.Photon(frequency=frequency)]:
            self.assertAlmostEqual(P.get_wavelength(), wavelength, 4)
            self.assertAlmostEqual(P.get_energy(), energy, 4)
            self.assertAlmostEqual(P.get_energy_eV(), energy_eV, 1)
            self.assertAlmostEqual(P.get_frequency()/1E15, frequency/1E15, 2)
