# ----------------------------------------------------------------------------------------------------- 
# PROPAGATOR: Scattering experiment simulator for spheres and customized object maps
# Please type 'help propagator()' for further information.
# -----------------------------------------------------------------------------------------------------
# Author:  Max Hantke - maxhantke@gmail.com
# -----------------------------------------------------------------------------------------------------
# All variables in SI units by default. Exceptions only if expressed by variable name.

import sys, ConfigParser, numpy, types, pickle, time, math, os

this_dir = os.path.dirname(os.path.realpath(__file__))

import logging
logger = logging.getLogger("Propagator")

# Initial configuration and importing propagator files
import config
config.init_configuration()
import imgutils,proptools
from source import Source
from sample import SampleMap,SampleSphere,SampleSpheroid
from detector import Detector

# Pythontools
from python_tools import gentools,cxitools,imgtools

class Input:
    """
    Input object which is the necessary argument of function 'propagator'
    It contains all the information about the experimental setup (objects 'source', 'sample', 'detector')

    """
    
    def __init__(self,configuration={}):
        """
        Function initializes input-object:
        ==================================
        Arguments:
        - configfile: Filename of configuration file. If not given variables are set to default values.

        """
        self.default_configuration = this_dir+"/data/default.conf"
        self.reconfigure(configuration)
        self._photon_changed = False
        self._detector_changed = False
    
    def reconfigure(self,configuration={}):
        """ 
        Function reconfigures Input subclasses based on the given configuration [self.configuration]
        """

        self.configuration = gentools.Configuration(configuration,self.default_configuration)

        C = self.configuration.confDict
        self.detector = Detector(parent=self,**C["detector"])
        self.source = Source(parent=self,**C["source"])

        if C["sample"]["sample_type"] == "uniform_sphere":
            self.sample = SampleSphere(parent=self,**C["sample"])
        elif C["sample"]["sample_type"] == "uniform_spheroid":
            self.sample = SampleSpheroid(parent=self,**C["sample"])
        elif C["sample"]["sample_type"] == "map3d":
            self.sample = SampleMap(parent=self,**C["sample"])
        else:
            logger.error("%s is not a valid sample type.")
            return

class Output:
    """
    OUTPUT of propagator provides user with results and functions for plotting.
    """
    def __init__(self,input_object):
        if not isinstance(input_object,Input):
            logger.error("Illegal input. Argument has to be of instance Input.")
            return
        
        self.input_object = input_object 
        logger.debug("Propagation started.")
        t_start = time.time()
        outdict = self.input_object.sample.propagate()
        self.amplitudes = outdict["amplitudes"]
        self.sample_euler_angle_0 = outdict.get("euler_angle_0",None)
        self.sample_euler_angle_1 = outdict.get("euler_angle_1",None)
        self.sample_euler_angle_2 = outdict.get("euler_angle_2",None)
        self.sample_diameter = outdict.get("sample_diameter",None)        
        t_stop = time.time()
        logger.debug("Propagation finished (time = %f sec)",t_stop-t_start)

    def get_intensity_pattern(self,i=0):
        """
        Returns 2-dimensional array with intensity values in photons per pixel (binned).
        """
        return self.input_object.detector.detect_photons(abs(self.amplitudes[i])**2)

    def get_real_space_image(self,i=0):
        A = self.amplitudes[i]
        A[numpy.isfinite(A)==False] = 0.
        return numpy.fft.fftshift(numpy.fft.ifftn(numpy.fft.fftshift(self.amplitudes[i])))

    def get_linear_oversampling_ratio(self):
        if self.input_object.sample.radius == None:
            return None
        else:
            pN = proptools.get_nyquist_pixel_size(self.input_object.detector.distance,self.input_object.source.photon.get_wavelength(),numpy.pi*self.input_object.sample.radius**2)
            pD = self.input_object.detector.get_pixel_size("binned")
            return pN/pD
            
    def get_full_period_edge_resolution(self):
        return proptools.get_max_crystallographic_resolution(self.input_object.source.photon.get_wavelength(),self.input_object.detector.get_minimum_center_edge_distance(),self.input_object.detector.distance)
