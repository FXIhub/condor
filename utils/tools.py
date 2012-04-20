import pylab, sys, numpy, types, pickle, time, math
#from constants import *
from propagator import *

def random_euler_angles():
    """
    Generates a triplet (phi, theta, psi) of random Euler angles.
    """
    r1,r2,r3 = pylab.random(3)
    q1 = pylab.sqrt(1.0-r1)*pylab.sin(2.0*pylab.pi*r2)
    q2 = pylab.sqrt(1.0-r1)*pylab.cos(2.0*pylab.pi*r2)
    q3 = pylab.sqrt(r1)*pylab.sin(2.0*pylab.pi*r3)
    q4 = pylab.sqrt(r1)*pylab.cos(2.0*pylab.pi*r3)
    e1 = math.atan2(2.0*(q1*q2+q3*q4), 1.0-2.0*(q2**2+q3**2))
    e2 = math.asin(2.0*(q1*q3-q4*q2))
    e3 = math.atan2(2.0*(q1*q4+q2*q3), 1.0-2.0*(q3**2+q4**2))
    return (e1,e2,e3)

def get_max_crystallographic_resolution(wavelength,min_detector_center_edge_distance,detector_distance):
    """
    Returns crystallographic resolution (full-period resolution at the closest edge)
    """
    return wavelength/pylab.sin(pylab.arctan(min_detector_center_edge_distance/detector_distance))
    
def get_nyquist_pixelsize(detector_distance,wavelength,particle_area):
    """
    Returns size of one Nyquist pixel on the detector in meter.
    """
    particle_radius = pylab.sqrt(particle_area/pylab.pi)
    return detector_distance * wavelength / (2*particle_radius)

def print_material_xray_properties(wavelength,thickness=1.0E-06,**margs):
    re = DICT_physical_constants['re']
    h = phy.DICT_physical_constants['h']
    c = phy.DICT_physical_constants['c']
    qe = phy.DICT_physical_constants['e']

    photon_energy_eV = h*c/wavelength/qe
    M = Material(photon_energy_eV,**margs)
    n = M.get_n(photon_energy_eV)
    print n
    f = M.get_f(photon_energy_eV)
    dn = n-1
    delta = -dn.real
    beta = -dn.imag
    phase_shift = 2*pylab.pi*thickness*delta/wavelength
    #mu_a = 2*re*wavelength*f2
    #s = 1/mu_a/n
    T = pylab.exp(-4*pylab.pi*beta/wavelength*thickness)
    print "BEAM:"
    print "Wavelength = %.2f nm ; Energy = %.0f eV" % (wavelength/1.0E-09,photon_energy_eV)
    print "SAMPLE DENSITY:"
    print "Mass denstity: %.3f mg/cm^3" % M.massdensity
    print "SAMPLE SCATTERING AND ABSORPTION PARAMETERS:"
    print "Scattering factor (real part) f1 = %e" % f.real
    print "Scattering factor (imag part) f2 = %e" % f.imag
    print "Refraction coefficient n = 1 - delta - i beta"
    print "delta = %f" % delta
    print "beta = %f" % beta
    print "Phaseshift / %.2f um = %f pi" % (thickness/1.0E-6,phase_shift/pylab.pi)
    #print "Atomic photoabsorption cross section: mu_a = %f re^2" % (mu_a/re**2)
    #print "Attenuation length (drop off to 1/e): s = %f um" % (s/1.0E-6)
    print "Transmission after %.2f um sample: T = %.1f percent " % (thickness/1.0E-06,T*100)
    #atomic photoabsorption cross section
