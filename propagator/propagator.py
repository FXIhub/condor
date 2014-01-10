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
import xcorepropagation,imgutils,proptools
import source,sample,detector

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
        self.configuration = gentools.Configuration(configuration,self.default_configuration)
        self.reconfigure()
        self._photon_changed = False
        self._detector_changed = False
    
    def reconfigure(self,configuration={}):
        """ 
        Function reconfigures Input subclasses based on the given configuration [self.configuration]
        """
        if configuration != {}:
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
        self.amplitudes = self.input_object.sample.propagate()
        t_stop = time.time()
        logger.debug("Propagation finished (time = %f sec)",t_stop-t_start)

    def get_intensity_pattern(self):
        """
        Returns 2-dimensional array with intensity values in photons per pixel (binned).
        """
        return abs(self.amplitudes)**2

    def get_real_space_image(self):
        A = self.amplitudes
        A[numpy.isfinite(A)==False] = 0.
        return numpy.fft.fftshift(numpy.fft.ifftn(numpy.fft.fftshift(self.amplitudes)))
