# -----------------------------------------------------------------------------------------------------
# CONDOR
# Simulator for diffractive single-particle imaging experiments with X-ray lasers
# http://xfel.icm.uu.se/condor/
# -----------------------------------------------------------------------------------------------------
# Copyright 2016 Max Hantke, Filipe R.N.C. Maia, Tomas Ekeberg
# Condor is distributed under the terms of the BSD 2-Clause License
# -----------------------------------------------------------------------------------------------------
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# 
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# -----------------------------------------------------------------------------------------------------
# General note:
# All variables are in SI units by default. Exceptions explicit by variable name.
# -----------------------------------------------------------------------------------------------------
import numpy
import nfft
import unittest

class TestNfft(unittest.TestCase):
    def setUp(self):
        self._size = 10
        self._decimals = 10
        self._coord_1d = numpy.linspace(-0.5, 0.5-1./self._size, self._size)
        self._2d_coord = [[3, 5], [8, 7]]
        self._3d_coord = [[3, 5, 2], [8, 7, 4]]
        self._4d_coord = [[3, 0, 2, 1], [3, 3, 1, 2]]

    def test_transformer_1d(self):
        a = numpy.random.random(self._size)
        t = nfft.Transformer(a, 20)
        ft_nfft = t.transform(self._coord_1d)
        ft_fftw = numpy.fft.fftshift(numpy.fft.fft(numpy.fft.fftshift(a)))
        numpy.testing.assert_almost_equal(ft_nfft, ft_fftw, decimal=self._decimals)

    def test_transformer_2d(self):
        a = numpy.random.random((self._size, )*2)
        t = nfft.Transformer(a, 5)
        ft_nfft = t.transform([[self._coord_1d[self._2d_coord[0][0]], self._coord_1d[self._2d_coord[0][1]]],
                               [self._coord_1d[self._2d_coord[1][0]], self._coord_1d[self._2d_coord[1][1]]]])
        ft_fftw = numpy.fft.fftshift(numpy.fft.fft2(numpy.fft.fftshift(a)))
        numpy.testing.assert_almost_equal(ft_nfft, ft_fftw[(self._2d_coord[0][0], self._2d_coord[1][0]), (self._2d_coord[0][1], self._2d_coord[1][1])])
        
    def test_nfft_1d(self):
        a = numpy.random.random(self._size)
        ft_nfft = nfft.nfft(a, self._coord_1d)
        ft_fftw = numpy.fft.fftshift(numpy.fft.fft(numpy.fft.fftshift(a)))
        numpy.testing.assert_almost_equal(ft_nfft, ft_fftw, decimal=self._decimals)

    def test_nfft_inplace_1d(self):
        a = numpy.random.random(self._size)
        ft_nfft = numpy.empty(self._size, dtype="complex128")
        nfft.nfft_inplace(a, self._coord_1d, ft_nfft)
        ft_fftw = numpy.fft.fftshift(numpy.fft.fft(numpy.fft.fftshift(a)))
        numpy.testing.assert_almost_equal(ft_nfft, ft_fftw, decimal=self._decimals)
            
    def test_nfft_2d(self):
        a = numpy.random.random((self._size, )*2)
        ft_nfft = nfft.nfft(a, [[self._coord_1d[self._2d_coord[0][0]], self._coord_1d[self._2d_coord[0][1]]],
                                [self._coord_1d[self._2d_coord[1][0]], self._coord_1d[self._2d_coord[1][1]]]])
        ft_fftw = numpy.fft.fftshift(numpy.fft.fft2(numpy.fft.fftshift(a)))
        numpy.testing.assert_almost_equal(ft_nfft, ft_fftw[(self._2d_coord[0][0], self._2d_coord[1][0]), (self._2d_coord[0][1], self._2d_coord[1][1])])

    def test_nfft_inplace_2d(self):
        a = numpy.random.random((self._size, )*2)
        ft_nfft = numpy.empty(2, dtype="complex128")
        nfft.nfft_inplace(a, [[self._coord_1d[self._2d_coord[0][0]], self._coord_1d[self._2d_coord[0][1]]],
                              [self._coord_1d[self._2d_coord[1][0]], self._coord_1d[self._2d_coord[1][1]]]], ft_nfft)
        ft_fftw = numpy.fft.fftshift(numpy.fft.fft2(numpy.fft.fftshift(a)))
        numpy.testing.assert_almost_equal(ft_nfft, ft_fftw[(self._2d_coord[0][0], self._2d_coord[1][0]), (self._2d_coord[0][1], self._2d_coord[1][1])])

    def test_nfft_3d(self):
        a = numpy.random.random((self._size, )*3)
        ft_nfft = numpy.empty(2, dtype="complex128")
        nfft.nfft_inplace(a, [[self._coord_1d[self._3d_coord[0][0]], self._coord_1d[self._3d_coord[0][1]], self._coord_1d[self._3d_coord[0][2]]],
                              [self._coord_1d[self._3d_coord[1][0]], self._coord_1d[self._3d_coord[1][1]], self._coord_1d[self._3d_coord[1][2]]]], ft_nfft)
        ft_fftw = numpy.fft.fftshift(numpy.fft.fftn(numpy.fft.fftshift(a)))
        numpy.testing.assert_almost_equal(ft_nfft, ft_fftw[(self._3d_coord[0][0], self._3d_coord[1][0]),
                                                           (self._3d_coord[0][1], self._3d_coord[1][1]),
                                                           (self._3d_coord[0][2], self._3d_coord[1][2])])

    def test_nfft_4d(self):
        size = 6
        coord_1d = numpy.linspace(-0.5, 0.5-1./size, size)
        a = numpy.random.random((size, )*4)
        ft_nfft = numpy.empty(2, dtype="complex128")
        nfft.nfft_inplace(a, [[coord_1d[self._4d_coord[0][0]], coord_1d[self._4d_coord[0][1]], coord_1d[self._4d_coord[0][2]], coord_1d[self._4d_coord[0][3]]],
                              [coord_1d[self._4d_coord[1][0]], coord_1d[self._4d_coord[1][1]], coord_1d[self._4d_coord[1][2]], coord_1d[self._4d_coord[1][3]]]],
                          ft_nfft)
        ft_fftw = numpy.fft.fftshift(numpy.fft.fftn(numpy.fft.fftshift(a)))
        numpy.testing.assert_almost_equal(ft_nfft, ft_fftw[(self._4d_coord[0][0], self._4d_coord[1][0]),
                                                           (self._4d_coord[0][1], self._4d_coord[1][1]),
                                                           (self._4d_coord[0][2], self._4d_coord[1][2]),
                                                           (self._4d_coord[0][3], self._4d_coord[1][3])])

    def test_transformer_direct(self):
        a = numpy.random.random(self._size)
        t = nfft.Transformer(a, 100)
        coordinates = -0.5+numpy.random.random(100)
        ft_nfft = t.transform(coordinates, use_direct=False)
        ft_nfft_direct = t.transform(coordinates, use_direct=True)
        numpy.testing.assert_almost_equal(ft_nfft, ft_nfft_direct, decimal=self._decimals)

    def test_nfft_direct(self):
        a = numpy.random.random(self._size)
        coordinates = -0.5+numpy.random.random(100)
        ft_nfft = nfft.nfft(a, coordinates, use_direct=False)
        ft_nfft_direct = nfft.nfft(a, coordinates, use_direct=True)
        numpy.testing.assert_almost_equal(ft_nfft, ft_nfft_direct, decimal=self._decimals)

    def test_nfft_inplace_direct(self):
        #a = numpy.random.random(self._size)
        a = numpy.random.random((4, )*4)
        #coordinates = -0.5+numpy.random.random(100)
        number_of_points = 100.
        coordinates = numpy.array(zip(-0.5 + numpy.random.random(number_of_points), -0.5 + numpy.random.random(number_of_points),
                                      -0.5 + numpy.random.random(number_of_points), -0.5 + numpy.random.random(number_of_points)))
        ft_nfft = numpy.empty(len(coordinates), dtype="complex128")
        nfft.nfft_inplace(a, coordinates, ft_nfft, use_direct=False)
        ft_nfft_direct = numpy.empty(len(coordinates), dtype="complex128")
        nfft.nfft_inplace(a, coordinates, ft_nfft_direct, use_direct=True)
        numpy.testing.assert_almost_equal(ft_nfft, ft_nfft_direct, decimal=3)

    def test_failures(self):
        a = numpy.random.random(self._size)
        self.assertRaises(TypeError, nfft.nfft, (a, "hej"))
        self.assertRaises(TypeError, nfft.nfft, ("hej", a))

if __name__ == "__main__":
    unittest.main()
    # my_tests = unittest.TestSuite()
    # my_tests.addTest(TestNfft("test_nfft_inplace_1d"))
    # my_tests.addTest(TestNfft("test_nfft_direct"))
    # my_tests.addTest(TestNfft("test_nfft_inplace_direct"))
    # unittest.TextTestRunner().run(my_tests)

