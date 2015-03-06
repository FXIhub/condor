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

import sys
sys.path.append("utils")
import numpy
import logging
logger = logging.getLogger("Condor")

import condortools
from variation import Variation
from spimage import PixelMask

# Pythontools
from python_tools import gentools,cxitools,imgtools

class Detector:
    """
    A subclass of the input object.
    Defines area detector.
    """
    def __init__(self,**kwargs):
        """
        Function initializes Detector object.
        =====================================
        Arguments:
        Keyword arguments (if not given variable is set to default value):

        - distance: Distance between detector and interaction point.
        - pixel_size: Edge length of square pixel (unbinned).
        - cx: Horizontal beam position in pixel. If argument is \'None\' or not given center is set to the middle. [None]
        - cy: Vertical beam position in pixel. If argument is \'None\' or not given center is set to the middle. [None]
        - center_variation: Variation of the center position. Either None, \'normal\', \'uniform\' or \'range\'. [None]
        - center_spread_x: If \'normal\', \'uniform\' or \'range\' center variation is specified spread defines width of the distribution in x.
        - center_spread_y: If \'normal\', \'uniform\' or \'range\' center variation is specified spread defines width of the distribution in y.
        - center_variation_n: I \'range\' center variation is specified this argument specifies the number of samples within the specified range.
        - binning: Number of binned pixels (binning x binning). [1]
        - noise: Noise that is put on the intensities when read. Either None, \'poisson\', \'normal\', \'uniform\' or \'normal_poisson\'. [None]
        - noise_spread: If \'normal\', \'uniform\' or \'normal_poisson\' noise is specified spread defines width of gaussian/uniform distribution.
        - parent: Input object that includes detector object. This variable is optional. [None]
        - saturation_level: Value at which detector pixels satutrate. [None]

        EITHER (default):
        - nx: Number of pixels in horizontal direction not including a potential gap.
        - ny: Number of pixels in vertical direction not including a potential gap.
        - x_gap_size_in_pixel: Width of central horizontal gap in pixel. [0]
        - y_gap_size_in_pixel: Width of central vertical gap in pixel. [0]
        - hole_diameter_in_pixel: Diameter of central hole in pixel. [0]
        
        OR:
        - mask_CXI_bitmask: bool wheter or not data is stored as a CXI bitmask. [False]
          + False: pixels with value 0 are invalid, pixels with value 1 are valid.
          + True: (pixels & spimage.PixelMask.PIXEL_IS_IN_MASK_DEFAULT) == 0 are the valid pixels
        
          EITHER
          - mask: 2d mask array. In the array 1 stands for valid pixel and 0 for an invalid pixel.
          OR:
          - mask_filename: path to HDF5 file where the mask is stored
          - mask_dataset: dataset path in file
        """

        # Check for required keyword arguments
        reqk = ["distance","pixel_size"]
        for k in reqk:
            if k not in kwargs.keys():
                logger.error("Cannot initialize Detector instance. %s is a necessary keyword." % k)
                return
        # Check for valid keyword arguments
        allk = ["parent",
                "distance","pixel_size",
                "x_gap_size_in_pixel","y_gap_size_in_pixel","hole_diameter_in_pixel",
                "noise","noise_spread","noise_variation_n",
                "center_variation","center_spread_x","center_spread_y","center_variation_n",
                "saturation_level","binning",
                "mask","mask_filename","mask_dataset",
                "cx","cy",
                "nx","ny"]
        self._unproc_kws = [k for k in kwargs.keys() if k not in allk]
        if len(self._unproc_kws) > 0:
            print self._unproc_kws
            logger.error("Detector object initialisation failed due to illegal keyword arguments.")
            return
        # Start initialisation            
        self._parent = kwargs.get("parent",None)
        self.distance = kwargs["distance"]
        self.pixel_size = kwargs["pixel_size"]
        gx = kwargs.get('x_gap_size_in_pixel',0)
        gy = kwargs.get('y_gap_size_in_pixel',0)
        hd = kwargs.get('hole_diameter_in_pixel',0)
        self.set_noise(noise=kwargs.get("noise",None),
                       noise_spread=kwargs.get("noise_spread",None),
                       noise_variation_n=kwargs.get("noise_variation_n",None))
        self.set_center_variation(center_variation=kwargs.get("center_variation",None),
                                  center_spread_x=kwargs.get("center_spread_x",None),
                                  center_spread_y=kwargs.get("center_spread_y",None),
                                  center_variation_n=kwargs.get("center_variation_n",None))
        self.saturation_level = kwargs.get('saturation_level',None)
        self._binning = kwargs.get('binning',1)
        if 'mask' in kwargs or ('mask_filename' in kwargs.keys() and "mask_dataset" in kwargs.keys()):
            if "mask" in kwargs:
                M = kwargs["mask"]
            else:
                f = h5py.File(kwargs["mask_filename"],"r")
                M = f["mask_dataset"][:]
                f.close()
            if kwargs.get("mask_CXI_bitmask",False): M = (M & PixelMask.PIXEL_IS_IN_MASK_DEFAULT) == 0
            self.init_mask(mask=M)
            self._cx_mean = kwargs.get('cx',(self._mask.shape[1]-1)/(2.*self._binning))
            self._cy_mean = kwargs.get('cy',(self._mask.shape[0]-1)/(2.*self._binning))              
            self._next()
        else:
            reqk = ["nx","ny"]
            for k in reqk:
                if k not in kwargs.keys():
                    logger.error("Cannot initialize Detector instance. %s is a necessary keyword if no mask is given." % k)
                    return
            self._nx = kwargs["nx"]
            self._ny = kwargs["ny"]
            self._cx_mean = kwargs.get('cx',(self._nx+gy-1)/2.)
            self._cy_mean = kwargs.get('cy',(self._ny+gx-1)/2.)
            self._next()
            self.init_mask(nx=kwargs["nx"],ny=kwargs["ny"],x_gap_size_in_pixel=gx,y_gap_size_in_pixel=gy,hole_diameter_in_pixel=hd,binning=self._binning)

    def set_noise(self, noise = None, noise_spread = None, noise_variation_n = None):
        self._noise = Variation(noise,noise_spread,noise_variation_n,number_of_dimensions=1,name="noise")

    def set_center_variation(self, center_variation = None, center_spread_x = None, center_spread_y = None, center_variation_n = None):
        self._center_variation = Variation(center_variation,[center_spread_y,center_spread_x],center_variation_n,number_of_dimensions=2,name="center position")

    def init_mask(self,**kwargs):
        """        
        Function initializes the detector mask.
        =======================================

        Arguments:
        
        Keyword arguments (if not given variable is set to default value):
        
        EITHER
        - mask: array that defines the mask (0: masked out, 1: not masked out)
        
        OR:
        - nx: horizontal dimension of the mask (inside the mask, without counting any gaps)
        - ny: vertical dimension of the mask (inside the mask, without counting any gaps)
        - x_gap_size_in_pixel: horizontal gap is generated with given width (in unit unbinned pixel)
        - y_gap_size_in_pixel: vertical gap is generated with given width (in unit unbinned pixel)
        - hole_diameter_in_pixel: holeis generated with given diameter (in unit unbinned pixel)

        """

        # init mask array
        if 'binning' in kwargs:
            self._binning = kwargs['binning']
        else:
            self._binning = 1
        if kwargs.get('mask',None) != None: 
            self._mask = kwargs['mask']
            self._nx = self._mask.shape[1]*self._binning
            self._ny = self._mask.shape[0]*self._binning
        elif 'nx' in kwargs and 'ny' in kwargs:
            self._nx = kwargs["nx"]
            self._ny = kwargs["ny"]
            Nx = kwargs['nx']/self._binning; Ny = kwargs['ny']/self._binning
            if 'x_gap_size_in_pixel' in kwargs:
                x_gap = kwargs['x_gap_size_in_pixel']/self._binning
                Ny += x_gap
            if 'y_gap_size_in_pixel' in kwargs: 
                y_gap = kwargs['y_gap_size_in_pixel']/self._binning
                Nx += y_gap
            self._mask = numpy.ones(shape=(Ny,Nx))
        else:
            logger.error("Either \'mask_array\' or \'nx\' and \'ny\' have to be specified.")
            return 
            
        # set pixels in gap to zero
        if 'x_gap_size_in_pixel' in kwargs:
            if kwargs['x_gap_size_in_pixel'] != 0:
                cy = numpy.ceil((self._mask.shape[0]-1)/2.)
                self._mask[cy-kwargs['x_gap_size_in_pixel']/2/self._binning:cy-kwargs['x_gap_size_in_pixel']/2/self._binning+kwargs['x_gap_size_in_pixel']/self._binning,:] = 0
        if 'y_gap_size_in_pixel' in kwargs:
            if kwargs['y_gap_size_in_pixel'] != 0:
                cx = numpy.ceil((self._mask.shape[1]-1)/2.)
                self._mask[:,cx-kwargs['y_gap_size_in_pixel']/2/self._binning:cx-kwargs['y_gap_size_in_pixel']/2/self._binning+kwargs['y_gap_size_in_pixel']/self._binning] = 0
        if 'hole_diameter_in_pixel' in kwargs:
            if kwargs['hole_diameter_in_pixel'] != 0:
                cy = (self._mask.shape[0]-1)/2.
                cx = (self._mask.shape[1]-1)/2.
                X,Y = numpy.meshgrid(numpy.arange(0,self._mask.shape[1],1.0),
                                     numpy.arange(0,self._mask.shape[0],1.0))
                X = X-cx
                Y = Y-cy
                R = numpy.sqrt(X**2 + Y**2)
                self._mask[R<=kwargs['hole_diameter_in_pixel']/(2.0*self._binning)] = 0
        

    def get_mask(self,intensities,output_bitmask=False):
        if output_bitmask:
            return self.get_bitmask(intensities)
        else:
            return numpy.array(self.get_bitmask(intensities) == 0,dtype="bool")
    
    def get_bitmask(self,intensities):
        M = numpy.zeros(shape=self._mask.shape,dtype="uint16")
        M[self._mask == 0] |= PixelMask.PIXEL_IS_MISSING
        if self.saturation_level is not None:
            M[intensities >= self.saturation_level] |= PixelMask.PIXEL_IS_SATURATED
        return M

    def set_cy(self,cy=None):
        """
        Function sets vertical center position:
        =======================================
        Arguments:
        cy : vertical center position in pixel. If argument is None or not given center is set to the middle.
        """
        if cy == None:
            self.cy = (self._mask.shape[0]-1)/2.
        else:
            self.cy = cy

    def set_cx(self,cx=None):
        """
        Function sets horicontal center position:
        =========================================
        Arguments:
        cx : horizontal center position in pixel. If argument is None or not given center is set to the middle.
        """
        if cx == None:
            self._cx_mean = (self._mask.shape[1]-1)/2.
        else:
            self._cx_mean = cx

    def _next(self):
        # Resolve relative center
        if self._cx_mean == "middle":
            cx = (self._nx-1)/2.
        else:
            cx = self._cx_mean
        if self._cy_mean == "middle":
            cy = (self._ny-1)/2.
        else:
            cy = self._cy_mean
        # Center variation
        self.cy,self.cx = self._center_variation.get([cy,cx])

    def get_cx(self,option="unbinned"):
        return self._get_c("x",option)

    def get_cy(self,option="unbinned"):
        return self._get_c("y",option)

    def _get_c(self,coord,option):
        # Choose coordinate
        if coord == "x":
            c = self.cx
        elif coord == "y":
            c = self.cy
        else:
            logger.error("No valid coordinate chosen.")
        # Take into account binning if required
        if option == 'unbinned':
            pass
        elif option == 'binned':
            if (c % 1) == 0.5:
                c = (c-(self._binning-1)/2.)/(1.0*self._binning)
            else:
                c = c/(1.0*self._binning)
        else:
            logger.error("No valid option chosen.")
        return c

    def get_minimum_center_edge_distance(self):
        cx = self.get_cx()
        icx = self._nx-cx
        cy = self.get_cy()
        icy = self._ny-cy
        return min([cx,icx,cy,icy])*self.get_pixel_size()

    def get_pixel_size(self,option='unbinned'):
        if option == 'unbinned':
            return self.pixel_size
        elif option == 'binned':
            return self.pixel_size*self._binning
        else:
            logger.error("No valid option chosen.")

    def get_pixel_solid_angle(self,option='unbinned'):
        return self.get_pixel_size(option)**2 / self.distance**2
    
    def generate_absqmap(self,**kwargs):
        X,Y = numpy.meshgrid(numpy.arange(self._mask.shape[1]),
                             numpy.arange(self._mask.shape[0]))
        # THIS CAST IS VERY IMPORTANT, in python A += B is not the same as A = A + B
        X = numpy.float64(X)
        Y = numpy.float64(Y)
        X -= self.get_cx('binned')
        Y -= self.get_cy('binned')
        p = self.get_pixel_size('binned')
        D = self.distance
        if "wavelength" in kwargs:
            w = kwargs["wavelength"]
        else:
            w = self._parent.source.photon.get_wavelength()
        return condortools.generate_absqmap(X,Y,p,D,w)

    def generate_qmap(self,**kwargs):
        X,Y = numpy.meshgrid(numpy.arange(self._mask.shape[1]),
                             numpy.arange(self._mask.shape[0]))
        # THIS CAST IS VERY IMPORTANT, in python A += B is not the same as A = A + B
        X = numpy.float64(X)
        Y = numpy.float64(Y)
        X -= self.get_cx('binned')
        Y -= self.get_cy('binned')
        p = self.get_pixel_size('binned')
        D = self.distance
        if "wavelength" in kwargs:
            w = kwargs["wavelength"]
        else:
            w = self._parent.source.photon.get_wavelength()
        if "euler_angle_0" in kwargs:
            E0 = kwargs["euler_angle_0"]
        else:
            E0 = self._parent.sample.euler_angle_0
        if "euler_angle_1" in kwargs:
            E1 = kwargs["euler_angle_1"]
        else:
            E1 = self._parent.sample.euler_angle_1
        if "euler_angle_2" in kwargs:
            E2 = kwargs["euler_angle_2"]
        else:
            E2 = self._parent.sample.euler_angle_2
        qmap = condortools.generate_qmap(X,Y,p,D,w,E0,E1,E2)
        nfft_scaled = kwargs.get("nfft_scaled",False)
        if nfft_scaled:
            qmap /= self.get_absq_max()/0.5*numpy.sqrt(2)
        return qmap
    
    def generate_qmap_ori(self,**kwargs):
        X,Y = numpy.meshgrid(numpy.arange(self._mask.shape[1]),
                             numpy.arange(self._mask.shape[0]))
        # THIS CAST IS VERY IMPORTANT, in python A += B is not the same as A = A + B
        X = numpy.float64(X)
        Y = numpy.float64(Y)
        X -= self.get_cx('binned')
        Y -= self.get_cy('binned')
        p = self.get_pixel_size('binned')
        D = self.distance
        if "wavelength" in kwargs:
            w = kwargs["wavelength"]
        else:
            w = self._parent.source.photon.get_wavelength()
        qmap = condortools.generate_qmap_ori(X,Y,p,D,w)
        qmap /= self.get_absq_max()/0.5*numpy.sqrt(2)
        return qmap

    def get_absqx_max(self,**kwargs):
        if "wavelength" in kwargs:
            w = kwargs["wavelength"]
        else:
            w = self._parent.source.photon.get_wavelength()
        x_max = max([self.get_cx('binned'),self._mask.shape[1]-1-self.get_cx('binned')]) * self.get_pixel_size('binned')
        R_Ewald = 2*numpy.pi/w
        phi = numpy.arctan2(x_max,self.distance)
        return 2 * R_Ewald * numpy.sin(phi/2.0)

    def get_absqy_max(self,**kwargs):
        if "wavelength" in kwargs:
            w = kwargs["wavelength"]
        else:
            w = self._parent.source.photon.get_wavelength()
        y_max = max([self.get_cy('binned'),self._mask.shape[0]-1-self.get_cy('binned')]) * self.get_pixel_size('binned')
        R_Ewald = 2*numpy.pi/w
        phi = numpy.arctan2(y_max,self.distance)
        return 2 * R_Ewald * numpy.sin(phi/2.0)

    def get_absqz_max(self,**kwargs):
        if "wavelength" in kwargs:
            w = kwargs["wavelength"]
        else:
            w = self._parent.source.photon.get_wavelength()
        absqx_max = self.get_absqx_max()
        absqy_max = self.get_absqy_max()
        w = self.source.photon.get_w()
        R_Ewald = 2*numpy.pi/w
        phi = numpy.arcsin(numpy.sqrt(absqx_max**2+absqy_max**2)/R_Ewald)
        return R_Ewald * (1-numpy.cos(phi))

    def get_absq_max(self,**kwargs):
        return max([self.get_absqx_max(**kwargs),self.get_absqy_max(**kwargs)])

    def get_real_space_resolution_element(self,**kwargs):
        dX = numpy.pi / self.get_absq_max(**kwargs)
        return dX

    def get_max_achievable_crystallographic_resolution(self,**kwargs):
        dx = 2 * numpy.pi / self.get_absqx_max(**kwargs)
        dy = 2 * numpy.pi / self.get_absqy_max(**kwargs)
        return [dx,dy]

    def detect_photons(self,I):
        I_det = self._noise.get(I)
        if self.saturation_level is not None:
            temp = I_det > self.saturation_level
            if temp.sum() > 0:
                I_det[temp] = self.saturation_level
        return I_det

