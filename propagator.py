# ----------------------------------------------------------------------------------------------------- 
# PROPAGATOR: Scattering experiment simulator for spheres and customized object maps
# Please type 'help propagator()' for further information.
# -----------------------------------------------------------------------------------------------------
# Author:  Max Hantke - maxhantke@gmail.com
# -----------------------------------------------------------------------------------------------------
# All variables in SI units by default. Exceptions only if expressed by variable name.

import pylab, sys, ConfigParser, numpy, types, pickle, time, math

import logging
logger = logging.getLogger("Propagator")

# For plotting
from matplotlib import rc
import matplotlib.pyplot as mpy
rc('text', usetex=True)
rc('font', family='serif')
mpy.rcParams['figure.figsize'] = 9, 9

# Initial configuration and importing propagator files
import config
config.init_configuration()
import xcorepropagation,imgutils,proptools
from source import *
from sample import *
from detector import *

# Pythontools
import gentools,cxitools,imgtools

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
        self.default_configuration = "conf/default.conf"
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

    #def get_nyquist_pixel_size(self):
    #    return proptools.get_nyquist_pixel_size(self.input_object.detector.distance,self.input_object.source.photon.get_wavelength(),self.input_object.sample.get_area())

    #def save_pattern_to_file(self,filename,**kwargs):
    #    """
    #    Function saves dataset to file of specified format.
    #    ===================================================
    #    
    #    Arguments:
    #    - filename: The file-format is specified using one of the following file-endings:
    #                - \'.h5\'
    #                - \'.png\'

        
    #    Keyword arguments:
    #    - log: True / False (default)
    #    - poisson: True / False (default)
    #    - colorscale: \'jet\' (default) / \'gray\'
    #    - use_spimage: True / False (default)

    #    """
    #    pattern = self.get_intensity_pattern()
    #    mask = self.input_object.detector.mask
    #    if 'poisson' in kwargs:
    #        if kwargs['poisson']:
    #            pattern = pylab.poisson(pattern)
    #    if 'log' in kwargs:
    #        if kwargs['log']:
    #            pattern = pylab.log10(pattern)
    #            pattern[pylab.isfinite(pattern)==False] = pattern[pylab.isfinite(pattern)].min()
    #    use_spimage = kwargs.get('use_spimage',False)
    #    colorscale = kwargs.get('colorscale','jet')
    #    if use_spimage:
    #        import spimage
    #        if filename[-3:]=='.h5':
    #            color = 0
    #        elif filename[-3:]=='.png':
    #            if colorscale  == 'gray':
    #                    color = 1
    #            elif colorscale  == 'jet':
    #                    color = 16
    #        else:
    #            logger.error("%s is not a valid fileformat for this function." % filename[-3:])
    #            return
    #        tmp_data = spimage.sp_image_alloc(pattern.shape[1],pattern.shape[0],1)
    #        tmp_data.image[:,:] = pattern[:,:]
    #        tmp_data.mask[:,:] = mask[:,:]
    #        spimage.sp_image_write(tmp_data,filename,0)
    #        spimage.sp_image_free(tmp_data)
    #    else:
    #        if filename[-4:]=='.png':
    #            pylab.imsave(filename,pattern*pylab.log10(mask*10),cmap=pylab.cm.get_cmap(colorscale))
    #        elif filename[-3:]=='.h5':
    #            import h5py
    #            f = h5py.File(filename,'w')
    #            pattern_ds = f.create_dataset('intensities', pattern.shape, pattern.dtype)
    #            pattern_ds[:,:] = pattern[:,:]
    #            amplitudes_ds = f.create_dataset('amplitudes', self.amplitudes.shape, self.amplitudes.dtype)
    #            amplitudes_ds[:,:] = self.amplitudes[:,:]
    #            f.close()
