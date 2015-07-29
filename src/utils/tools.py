# -----------------------------------------------------------------------------------------------------
# CONDOR
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

import numpy, sys, numpy, types, pickle, time, math
import icosahedron
 
import logging
logger = logging.getLogger("Condor")
from log import log

def print_material_xray_properties(wavelength,thickness=1.0E-06,**margs):
    #r_0 = constants.value("classical electron radius")
    h =  constants.h
    c =  constants.c
    qe = constants.e

    photon_energy_eV = h*c/wavelength/qe
    M = Material(photon_energy_eV,**margs)
    n = M.get_n(photon_energy_eV)
    print n
    f = M.get_f(photon_energy_eV)
    dn = n-1
    delta = -dn.real
    beta = -dn.imag
    phase_shift = 2*numpy.pi*thickness*delta/wavelength
    #mu_a = 2*re*wavelength*f2
    #s = 1/mu_a/n
    T = numpy.exp(-4*numpy.pi*beta/wavelength*thickness)
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
    print "Phaseshift / %.2f um = %f pi" % (thickness/1.0E-6,phase_shift/numpy.pi)
    #print "Atomic photoabsorption cross section: mu_a = %f re^2" % (mu_a/re**2)
    #print "Attenuation length (drop off to 1/e): s = %f um" % (s/1.0E-6)
    print "Transmission after %.2f um sample: T = %.1f percent " % (thickness/1.0E-06,T*100)
    #atomic photoabsorption cross section



class Input:
    """
    The Input object that holds all necessary information for the simulation experiment. After initialization the configuration is saved to the variable configuration.confDict.

    :param configuration: Either a dictionary or the location of the configuration file. Missing but necessary arguments will be set to default values as specified in *default.conf*.
    
    """
    
    def __init__(self,configuration=None):
        C_raw = config.Configuration(configuration).confDict
        self.confDict = {}

        # Source
        self.confDict["source"] = config.Configuration(C_raw, {"source": get_default_source_conf()}).confDict["source"]
        self.source = Source(**self.confDict["source"])

        # Sample
        self.confDict["sample"] = config.Configuration(C_raw, {"sample": get_default_sample_conf()}).confDict["sample"]
        self.sample = Sample(**self.confDict["sample"])

        # Particles
        for k in [k for k in C_raw.keys() if k.startswith("particle")]:
            if k.endswith("sphere"):
                self.confDict[k] = config.Configuration({"particle":C_raw[k]}, {"particle": get_default_particle_uniform_sphere_conf()}).confDict["particle"]
                P = ParticleModelSphere(**self.confDict[k])
            elif k.endswith("spheroid"):
                self.confDict[k] = config.Configuration({"particle":C_raw[k]}, {"particle": get_default_particle_uniform_spheroid_conf()}).confDict["particle"]
                P = ParticleModelSpheroid(**self.confDict[k])
            elif k.endswith("map"):
                self.confDict[k] = config.Configuration({"particle":C_raw[k]}, {"particle": get_default_particle_map_conf()}).confDict["particle"]
                P = ParticleModelMap(**self.confDict[k])
            elif k.endswith("molecule"):
                self.confDict[k] = config.Configuration({"particle":C_raw[k]}, {"particle": get_default_particle_molecule_conf()}).confDict["particle"]
                P = ParticleModelMolecule(**self.confDict[k])
            else:
                log(logger.error,"Particle model for %s is not implemented." % k)
                sys.exit(1)
            self.sample.particle_models.append(P)
            
        # Detector
        self.confDict["detector"] = config.Configuration(C_raw, {"detector": get_default_detector_conf()}).confDict["detector"]
        self.detector = Detector(**self.confDict["detector"])

    def _write_conf(self, filename):
        log(logger.debug,"Writing configuration to: %s" % filename) 
        C = config.Configuration(self.confDict)
        C.write_to_file(filename)

class Output:
    """
    An instance of the Output object is initialized with an instance of the Input object and initiates the simulation of the diffraction data.
    After completion the instance holds the results and methods to access and interpret them.

    """
    def __init__(self,input):
        if not isinstance(input,Input):
            log(logger.error,"Illegal input. Argument has to be of instance Input.")
            sys.exit(1)
        self.input_object = input 
        self.experiment = Experiment(input.source,input.sample,input.detector)
        log(logger.debug,"Propagation started.")
        t_start = time.time()
        self.outdict = self.experiment.propagate()
        # General variables
        self.fourier_pattern   = self.outdict["fourier_pattern"]
        self.intensity_pattern = self.outdict["intensity_pattern"]
        self.mask              = self.outdict["mask_binary"]
        self.bitmask           = self.outdict["mask"]
        self.N = len(self.fourier_pattern)
        t_stop = time.time()
        log(logger.debug,"Propagation finished (time = %f sec)" % (t_stop-t_start))
        confout = "./condor.confout"
        self.input_object._write_conf(confout)
                
    def get_mask(self,i=0,output_bitmask=False):
        """
        Returns 2-dimensional array with bit mask values (binned).

        :param i: Index of the image that you want to obtain.
        :param output_bitmask: If True mask is written as CXI bit mask (16 bits, encoding see below).

        Bit code (16 bit integers) from spimage.PixelMask:
        PIXEL_IS_PERFECT = 0
        PIXEL_IS_INVALID = 1
        PIXEL_IS_SATURATED = 2
        PIXEL_IS_HOT = 4
        PIXEL_IS_DEAD = 8
        PIXEL_IS_SHADOWED = 16
        PIXEL_IS_IN_PEAKMASK = 32
        PIXEL_IS_TO_BE_IGNORED = 64
        PIXEL_IS_BAD = 128
        PIXEL_IS_OUT_OF_RESOLUTION_LIMITS = 256
        PIXEL_IS_MISSING = 512
        PIXEL_IS_NOISY = 1024
        PIXEL_IS_ARTIFACT_CORRECTED = 2048
        PIXEL_FAILED_ARTIFACT_CORRECTION = 4096
        PIXEL_IS_PEAK_FOR_HITFINDER = 8192
        PIXEL_IS_PHOTON_BACKGROUND_CORRECTED = 16384

        """
        return self.input_object.detector.get_mask(self.get_intensity_pattern(i),output_bitmask)

    def get_real_space_image(self,i=0):
        """
        Returns 2-dimensional array of back-propagated real space image from the diffraction fourier pattern.

        :param i: Index of the image that you want to obtain.

        """       
        A = self.fourier_pattern[i]
        A[numpy.isfinite(A)==False] = 0.
        return numpy.fft.fftshift(numpy.fft.ifftn(numpy.fft.fftshift(self.fourier_pattern[i])))

    def get_linear_sampling_ratio(self):
        """
        Returns the linear sampling ratio :math:`o` of the diffraction pattern:

        | :math:`o=\\frac{D\\lambda}{dp}` 

        | :math:`D`: Detector distance
        | :math:`p`: Detector pixel size (edge length)
        | :math:`\\lambda`: Photon wavelength 
        | :math:`d`: Sample diameter

        """       
        
        if self.input_object.sample.radius == None:
            return None
        else:
            pN = utils.diffraction.get_nyquist_pixel_size(self.input_object.detector.distance,self.input_object.source.photon.get_wavelength(),numpy.pi*self.input_object.sample.radius**2)
            pD = self.input_object.detector.get_pixel_size("binned")
            return pN/pD

    def write(self,filename="out.cxi",output="all"):
        if filename[-len(".cxi"):] != ".cxi":
            log(logger.error,"Illegal file format chosen.")
            return
        allout = output == "all"
        W = utils.log.CXIWriter(filename,self.N)

        def add(d0,d1,i):
            for k,v in d0.items():
                if k[0] == "_":
                    pass
                elif isinstance(v,dict):
                    d1[k] = {}
                    add(d0[k],d1[k],i)
                else:
                    d1[k] = d0[k][i]

        for i in range(self.N):
            O = {}
            if allout or "intensity_pattern" in output:
                O["intensity_pattern"] = self.intensity_pattern[i][:,:]
            if allout or "fourier_space_image" in output:
                O["fourier_pattern"] = self.fourier_pattern[i][:,:]
            if allout or "real_space_image" in output:
                O["real_space_image"] = self.get_real_space_image(i)
            add(self.outdict,O,i)
            W.write(O,i=i)
        W.close()


def get_default_source_conf():
    return config.Configuration(this_dir+"/data/source.conf").confDict["source"]    

def get_default_sample_conf():
    return config.Configuration(this_dir+"/data/sample.conf").confDict["sample"]    

def get_default_detector_conf():
    return config.Configuration(this_dir+"/data/detector.conf").confDict["detector"]    

def get_default_particle_uniform_sphere_conf():
    return config.Configuration(this_dir+"/data/particle_uniform_sphere.conf").confDict["particle"]    

def get_default_particle_uniform_spheroid_conf():
    return config.Configuration(this_dir+"/data/particle_uniform_spheroid.conf").confDict["particle"]    

def get_default_particle_map_conf():
    return config.Configuration(this_dir+"/data/particle_map.conf").confDict["particle"]    

def get_default_particle_molecule_conf():
    return config.Configuration(this_dir+"/data/particle_molecule.conf").confDict["particle"]    
