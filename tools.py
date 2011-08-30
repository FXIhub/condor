# Import library packages
#------------
import pylab, sys, numpy, types, pickle, time, math


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

def get_material_xray_properties(wavelength,thickness=1.0E-06,**margs):
    re = DICT_physical_constants['re']
    inp = Input()
    inp.source.wavelength = wavelength
    inp.set_sample_homogeneous_sphere(1,**margs)
    [nf1,nf2,n] = inp.sample.material.get_f_times_n0()
    f1 = nf1/n
    f2 = nf2/n
    n1 = 1-re/2/pylab.pi*wavelength**2*nf1
    n2 = -re/2/pylab.pi*wavelength**2*nf2
    delta = -n1+1
    beta = -n2
    phase_shift = 2*pylab.pi*thickness*delta/wavelength
    mu_a = 2*re*wavelength*f2
    s = 1/mu_a/n
    T = pylab.exp(-mu_a*n*thickness)
    print "BEAM:"
    print "Wavelength = %.2f nm ; Energy = %.0f eV" % (wavelength/1.0E-09,convert_photon(wavelength,"m","eV"))
    print "SAMPLE DENSITY:"
    print "Mass denstity: %.3f mg/cm^3" % inp.sample.material.massdensity
    print "Average atom denstity: %.3f N_A/cm^3" % (n/1.0E6/1.6E23)
    print "SAMPLE SCATTERING AND ABSORPTION PARAMETERS:"
    print "Scattering factor (real part) f1 = %f" % f1
    print "Scattering factor (imag part) f2 = %f" % f2
    print "Refraction coefficient n = 1 - delta - i beta = %f - i %f" % (n1,-n2)
    print "delta = %f" % delta
    print "beta = %f" % beta
    print "Phaseshift / %.2f mu = %f pi" % (thickness/1.0E-6,phase_shift/pylab.pi)
    print "Atomic photoabsorption cross section: mu_a = %f re^2" % (mu_a/re**2)
    print "Attenuation length (drop off to 1/e): s = %f mum" % (s/1.0E-6)
    print "Transmission after %.2f mu sample: T = %.1f percent " % (thickness/1.0E-06,T*100)
    #atomic photoabsorption cross section
