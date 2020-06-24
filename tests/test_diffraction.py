import unittest
from condor.utils import diffraction

class TestCaseDiffraction(unittest.TestCase):
    def test_recolution(self):
        # Seibert et al. 2011, Single mimivirus particles intercepted and imaged with an X-ray laser
        detector_distance = 0.564
        pixel_center_distance = 75E-6 * 512
        wavelength = 0.69E-9
        res_expected = 10.2E-9
        resel_expected = res_expected / 2.
        res = diffraction.crystallographic_resolution(wavelength, pixel_center_distance, detector_distance)
        resel = diffraction.resolution_element(wavelength, pixel_center_distance, detector_distance)
        self.assertAlmostEqual(res/1E-9 , res_expected/1E-9, 1)
        self.assertAlmostEqual(resel/1E-9 , resel_expected/1E-9, 1)

    def test_nyquist_pixel(self):
        # wavelength / particle_size = speckle_size / detector_distance
        wavelength = 1E-9
        detector_distance = 1.
        particle_size = 100E-9
        nypx_expected  = 0.01
        nypx = diffraction.nyquist_pixel_size(wavelength, detector_distance, particle_size)
        self.assertAlmostEqual(nypx/1E-3 , nypx_expected/1E-3, 1)
