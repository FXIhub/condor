#!/usr/bin/env python
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

import sys, numpy, types, time, os

this_dir = os.path.dirname(os.path.realpath(__file__))

sys.path.insert(0,os.path.join(this_dir, "utils"))

import logging
logger = logging.getLogger("Condor")
import utils
import utils.log
from utils.log import log 

# Initial configuration and importing Condor modules
import config
config.init_configuration()
from source import Source
from sample import Sample
from detector import Detector
from propagator import Propagator 
from particle_sphere import ParticleSpeciesSphere
from particle_spheroid import ParticleSpeciesSpheroid
from particle_map import ParticleSpeciesMap
from particle_molecule import ParticleSpeciesMolecule


class Input:
    """
    The Input object that holds all necessary information for the simulation experiment. After initialization the configuration is saved to the variable configuration.confDict.

    :param configuration: Either a dictionary or the location of the configuration file. Missing but necessary arguments will be set to default values as specified in *default.conf*.
    
    """
    
    def __init__(self,configuration={}):
        default_configuration = this_dir+"/data/default.conf"
        C = config.Configuration(configuration,default_configuration).confDict
        self.source = Source(**C["source"])
        self.sample = Sample(**C["sample"])
        for k in [k for k in C.keys() if "particle" in k]:
            if "particle_species" not in C[k]:
                log(logger.error,"No particle species defined for %s" % k)
                sys.exit(1)
            t = C[k]["particle_species"]
            if t == "uniform_sphere":
                P = ParticleSpeciesSphere(**C[k])
            elif t == "uniform_spheroid":
                P = ParticleSpeciesSpheroid(**C[k])
            elif t == "map3d":
                P = ParticleSpeciesMap(**C[k])
            elif t == "molecule":
                P = ParticleSpeciesMolecule(**C[k])
            else:
                log(logger.error,"ParticleSpecies class for particle_species=%s is not implemented." % t)
                sys.exit(1)
            self.sample.particle_species.append(P)        
        self.detector = Detector(**C["detector"])     

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
        self.propagator = Propagator(input.source,input.sample,input.detector)
        log(logger.debug,"Propagation started.")
        t_start = time.time()
        self.outdict = self.propagator.propagate()
        # General variables
        self.fourier_pattern   = self.outdict["fourier_pattern"]
        self.intensity_pattern = self.outdict["intensity_pattern"]
        self.mask              = self.outdict["mask_binary"]
        self.bitmask           = self.outdict["mask"]
        self.N = len(self.fourier_pattern)
        t_stop = time.time()
        log(logger.debug,"Propagation finished (time = %f sec)" % (t_stop-t_start))
    
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
            
    def get_full_period_edge_resolution(self):
        """
        Returns the full-period resolution :math:`R` at the edge of the detector in meter.

        | :math:`R=\\lambda / \\sin(\\arctan(Y/D))`

        | :math:`\\lambda`: Photon wavelength
        | :math:`Y`: Minimum distance between the beam axis and an edge of the detector.
        | :math:`D`: Detector distance

        """
        return utils.diffraction.get_max_crystallographic_resolution(self.input_object.source.photon.get_wavelength(),self.input_object.detector.get_minimum_center_edge_distance(),self.input_object.detector.distance)

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


