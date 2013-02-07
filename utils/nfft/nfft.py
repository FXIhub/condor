"""
Module that performs some nonequispaced ffts using the nfft C library.

This module is made for simulation of diffraction
patterns and is not intended for general use.
"""

import pylab
import nfft_c as _nfft_c

def nfft1d(coordinates, vin):
    """
    Performs a 1d nonequispaced fft.

    Input:
    coordinates : array of coordinates where the fft is evaluated. Withing the range [-0.5,0.5]
    vin         : input 1d array to transform

    Output:
    vout        : values of the transform at the positions specified in coordinates
    """
    vout = pylab.zeros(len(coordinates), dtype='complex128')

    vin.dtype = 'float64'
    vout.dtype = 'float64'

    _nfft_c.nfft(coordinates, vin, vout)

    vin.dtype = 'complex128'
    vout.dtype = 'complex128'
    
    return vout

def nfft3d(coordinates, vin):
    """
    Performs a 3d nonequispaced fft.

    Input:
    coordinates : array of triplets of coordinates where the fft is evaluated. Within the range [-0.5,0.5].
    vin         : input 3d array to transform

    Output:
    vout : values of the transform at the positions specified in coordinates
    """

    good_values_1d = (coordinates.flatten() <= 0.5) & (coordinates.flatten() >= -0.5)
    if sum(good_values_1d)!=len(good_values_1d):
        print "ERROR: nfft coordinates have to lie between -0.5 and 0.5."
        print 'number of good values = %d (%d)' % (sum(good_values_1d), len(good_values_1d))
        return []
    if vin.dtype != "complex128":
        vin = pylab.complex128(vin)
        print "Convert to complex128"
    vout_real = pylab.zeros(len(coordinates)/3, dtype='float64')
    vout_imag = pylab.zeros(len(coordinates)/3, dtype='float64')
    vin_real = pylab.array(vin.real, dtype='float64')
    vin_imag = pylab.array(vin.imag, dtype='float64')
    _nfft_c.nfft3(coordinates, vin_real, vin_imag, vout_real, vout_imag)
    vout = pylab.zeros(len(coordinates)/3, dtype='complex128')
    vout.real[:] = vout_real[:]
    vout.imag[:] = vout_imag[:]
    #_nfft_c.nfft3_finalize(coordinates, vin_real, vin_imag, vout_real, vout_imag)
    return vout

def test_nfft3d(N_sample=1000,N_coordinates=1000):
    cube = pylab.ones(shape=(N_sample,N_sample,N_sample))
    #print "2 %f" % (psutil.avail_phymem()/(1.0*mem0))
    x = pylab.arange(0,N_coordinates,1)/(1.0*(N_coordinates-1))-0.5
    #print "3 %f" % (psutil.avail_phymem()/(1.0*mem0))
    X,Y = pylab.meshgrid(x,x)
    #print "4 %f" % (psutil.avail_phymem()/(1.0*mem0))
    #Z = pylab.zeros(shape=(N,N))
    print N_coordinates
    print len(X.flatten())
    coordinates = pylab.zeros(shape=(N_coordinates*N_coordinates,3))
    #print "5 %f" % (psutil.avail_phymem()/(1.0*mem0))
    coordinates[:,1] = Y.flatten()
    #print "6 %f" % (psutil.avail_phymem()/(1.0*mem0))
    coordinates[:,2] = X.flatten()
    #print "7 %f" % (psutil.avail_phymem()/(1.0*mem0))
    coordinates = coordinates.flatten()
    #print "8 %f" % (psutil.avail_phymem()/(1.0*mem0))
    return nfft3d(coordinates,cube)

