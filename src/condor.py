# ----------------------------------------------------------------------------------------------------- 
# CONDOR 
# Simulator for diffractive single-particle imaging experiments with X-ray lasers
# http://xfel.icm.uu.se/condor/
# ----------------------------------------------------------------------------------------------------- 
# Copyright 2014 Max Hantke, Filipe R.N.C. Maia, Tomas Ekeberg
# Condor is distributed under the terms of the GNU General Public License
# ----------------------------------------------------------------------------------------------------- 
# General note:
#  All variables are in SI units by default. Exceptions explicit by variable name.
# ----------------------------------------------------------------------------------------------------- 

import sys, ConfigParser, numpy, types, pickle, time, math, os

this_dir = os.path.dirname(os.path.realpath(__file__))

sys.path.append(os.path.join(this_dir, "utils/python_tools"))

import logging
logger = logging.getLogger("Condor")

# Initial configuration and importing Condor files
import config
config.init_configuration()
import imgutils,condortools
from python_tools import cxitools
from source import Source
from sample_multiple import SampleMultiple
from sample_sphere import SampleSphere
from sample_spheroid import SampleSpheroid
from sample_map import SampleMap
from detector import Detector

# Pythontools
from python_tools import gentools,cxitools,imgtools

class Input:
    """
    The Input object that holds all necessary information for the experiment that shall be simulated. After initialization the configuration is saved to the variable configuration.confDict.

    :param configuration: Either a dictionary or the location of the configuration file. Missing but necessary arguments will be set to default values as specified in *default.conf*.
    
    """
    
    def __init__(self,configuration={}):
        self.default_configuration = this_dir+"/data/default.conf"
        self._reconfigure(configuration)
        self._photon_changed = False
        self._detector_changed = False
    
    def _reconfigure(self,configuration={}):
        self.configuration = gentools.Configuration(configuration,self.default_configuration)

        C = self.configuration.confDict
        self.detector = Detector(parent=self,**C["detector"])
        self.source = Source(parent=self,**C["source"])
        
        samples = [s for s in C.keys() if "sample" in s]
        self.sample = SampleMultiple(parent=self,**C[samples[0]])
        for s in samples:
            Cs = C[s]
            if Cs["sample_type"] == "uniform_sphere":
                S = SampleSphere(parent=self,**Cs)
            elif Cs["sample_type"] == "uniform_spheroid":
                S = SampleSpheroid(parent=self,**Cs)
            elif Cs["sample_type"] == "map3d":
                S = SampleMap(parent=self,**Cs)
            else:
                logger.error("%s is not a valid sample type.")
                return
            self.sample.add(S,Cs.get("translation_vector",None))

class Output:
    """
    An instance of the Output object is initialized with an instance of the Input object and initiates the simulation of the diffraction data.
    After completion the instance holds the results and methods to access and interpret them.

    """
    def __init__(self,input):
        if not isinstance(input,Input):
            logger.error("Illegal input. Argument has to be of instance Input.")
            return
        
        self.input_object = input 
        logger.debug("Propagation started.")
        t_start = time.time()
        outdict = self.input_object.sample.propagate()
        # General variables
        self.amplitudes = outdict["amplitudes"]
        self.N = len(self.amplitudes)
        self.sample_diameter = outdict.get("sample_diameter",None)
        self.cx = outdict.get("cx")
        self.cy = outdict.get("cy")
        self.cxXxX = outdict.get("cxXxX")
        self.cyXxX = outdict.get("cyXxX")
        self.beam_intensity = outdict["beam_intensity"]
        self._optional_outputs = ["F0","dX3","grid","qmap3d",
                                  "sample_spheroid_diameter_a","sample_spheroid_diameter_c","sample_spheroid_theta","sample_spheroid_phi","sample_spheroid_flattening","euler_angle_0","euler_angle_1","euler_angle_2"]
        for s in self._optional_outputs:
            if s in outdict:
                exec "self.%s = outdict[\"%s\"]" % (s,s)
        self._intensities = {}
        t_stop = time.time()
        logger.debug("Propagation finished (time = %f sec)",t_stop-t_start)
    
    def get_intensity_pattern(self,i=0):
        """
        Returns 2-dimensional array with intensity values in the unit photons per pixel (binned).

        :param i: Index of the image that you want to obtain.

        """
        if i not in self._intensities:
            self._intensities[i] = self.input_object.detector.detect_photons(abs(self.amplitudes[i])**2)
        return self._intensities[i]

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
        Returns 2-dimensional array of back-propagated real space image from the diffraction amplitudes.

        :param i: Index of the image that you want to obtain.

        """       
        A = self.amplitudes[i]
        A[numpy.isfinite(A)==False] = 0.
        return numpy.fft.fftshift(numpy.fft.ifftn(numpy.fft.fftshift(self.amplitudes[i])))

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
            pN = condortools.get_nyquist_pixel_size(self.input_object.detector.distance,self.input_object.source.photon.get_wavelength(),numpy.pi*self.input_object.sample.radius**2)
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
        return condortools.get_max_crystallographic_resolution(self.input_object.source.photon.get_wavelength(),self.input_object.detector.get_minimum_center_edge_distance(),self.input_object.detector.distance)

    def write(self,filename="out.cxi",output=["intensity_pattern","fourier_space_image"]):
        if filename[-len(".cxi"):] != ".cxi":
            logger.error("Illegal file format chosen.")    
            return
        allout = output == "all"
        W = cxitools.CXIWriter(filename,self.N,logger)
        for i in range(0,self.N):
            O = {}
            if allout or "intensity_pattern" in output:
                O["/entry_1/data_1/data"] = self.get_intensity_pattern(i)
            if allout or "real_space_image" in output:
                O["/entry_1/data_1/real_space_image"] = self.get_real_space_image(i)
            if allout or "fourier_space_image" in output:
                O["/entry_1/data_1/fourier_space_image"] = self.amplitudes[i][:,:]
            if allout or "bitmask_image" in output:
                O["/entry_1/data_1/mask"] = self.get_mask(i,output_bitmask=True)
            if allout or "mask_image" in output:
                O["/entry_1/data_1/mask_binary"] = self.get_mask(i,output_bitmask=False)
            if allout or "sample_diameter" in output:
                O["/entry_1/data_1/sample_diameter"] = self.sample_diameter[i]
            if allout or "beam_intensity" in output:
                O["/entry_1/data_1/beam_intensity"] = self.beam_intensity[i]
            if allout or "intensity_pattern_center" in output:
                O["entry_1/data_1/intensity_pattern_center_x"] = self.cx[i] 
                O["entry_1/data_1/intensity_pattern_center_y"] = self.cy[i]
                O["entry_1/data_1/intensity_pattern_center_x_binned"] = self.cxXxX[i] 
                O["entry_1/data_1/intensity_pattern_center_y_binned"] = self.cyXxX[i]
            for s in self._optional_outputs:
                if allout or s in output:
                    if not hasattr(self,s):
                        if not allout:
                            logger.warning("Cannot find output: %s",s)
                    else:
                        exec ("O[\"entry_1/data_1/%s\"] = self.%s[i]" % (s,s))
            W.write(O,i=i)
        W.close()
        
        
