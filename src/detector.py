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
        allk = ["distance","pixel_size",
                "x_gap_size_in_pixel","y_gap_size_in_pixel","hole_diameter_in_pixel",
                "noise","noise_spread","noise_variation_n",
                "center_variation","center_spread_x","center_spread_y","center_variation_n","center_spread_limit",
                "saturation_level","binning",
                "mask","mask_filename","mask_dataset",
                "cx","cy",
                "nx","ny"]
        self._unproc_kws = [k for k in kwargs.keys() if k not in allk]
        if len(self._unproc_kws) > 0:
            print self._unproc_kws
            logger.error("Detector object initialisation failed due to illegal keyword arguments.")
            exit(1)
        # Start initialisation            
        self.distance = kwargs["distance"]
        self.pixel_size = kwargs["pixel_size"]
        gx = kwargs["x_gap_size_in_pixel"]
        gy = kwargs["y_gap_size_in_pixel"]
        hd = kwargs["hole_diameter_in_pixel"]
        self.set_noise(noise=kwargs["noise"],
                       noise_spread=kwargs.get("noise_spread",None),
                       noise_variation_n=kwargs.get("noise_variation_n",None))
        self.set_center_variation(center_variation=kwargs["center_variation"],
                                  center_spread_x=kwargs["center_spread_x"],
                                  center_spread_y=kwargs["center_spread_y"],
                                  center_variation_n=kwargs.get("center_variation_n",None),
                                  center_spread_limit=kwargs["center_spread_limit"])
        self.saturation_level = kwargs["saturation_level"]
        self._binning = kwargs["binning"]
        if "mask" in kwargs or  ("mask_filename" in kwargs and "mask_dataset" in kwargs):
            if "mask" in kwargs:
                M = kwargs["mask"]
            else:
                with h5py.File(kwargs["mask_filename"],"r") as f:
                    M = f["mask_dataset"][:]
            if kwargs.get("mask_CXI_bitmask",False): M = (M & PixelMask.PIXEL_IS_IN_MASK_DEFAULT) == 0
            self.init_mask(mask=M)
            self._cx_mean = kwargs.get('cx',(self._mask.shape[1]-1)/(2.*self._binning))
            self._cy_mean = kwargs.get('cy',(self._mask.shape[0]-1)/(2.*self._binning))              
        else:
            reqk = ["nx","ny"]
            for k in reqk:
                if k not in kwargs.keys():
                    logger.error("Cannot initialize Detector instance. %s is a necessary keyword if no mask is given." % k)
                    exit(1)
            self._nx = kwargs["nx"]
            self._ny = kwargs["ny"]
            self._cx_mean = kwargs.get('cx',(self._nx+gy-1)/2.)
            self._cy_mean = kwargs.get('cy',(self._ny+gx-1)/2.)
            self.init_mask(nx=kwargs["nx"],ny=kwargs["ny"],x_gap_size_in_pixel=gx,y_gap_size_in_pixel=gy,hole_diameter_in_pixel=hd,binning=self._binning)

    def set_noise(self, noise, noise_spread, noise_variation_n):
        self._noise = Variation(noise,noise_spread,noise_variation_n,number_of_dimensions=1,name="noise")

    def set_center_variation(self, center_variation, center_spread_x, center_spread_y, center_variation_n, center_spread_limit):
        self._center_variation = Variation(center_variation,[center_spread_y,center_spread_x],center_variation_n,number_of_dimensions=2,name="center position")
        self._center_spread_limit = center_spread_limit
        
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

    def get_nx_binned(self):
        return self._mask.shape[1]
        
    def get_ny_binned(self):
        return self._mask.shape[0]

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

    def get_next(self):
        O = {}
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
        ready = False 
        while not ready:
            cy_now, cx_now = self._center_variation.get([cy,cx])
            ready = True
            if self._center_spread_limit > 0:
                if (self._center_spread_limit < abs(cx_now - self._cx_mean)*2) and (self._center_spread_limit < abs(cy_now - self._cy_mean)*2):
                    ready = False
        O["cx_unbinned"] = cx_now
        O["cy_unbinned"] = cy_now
        O["cx_binned"] = self._get_c(cx_now,"binned")
        O["cy_binned"] = self._get_c(cy_now,"binned")
        O["solid_angle_binned_pixel"] = self.get_pixel_solid_angle("binned")
        O["nx_binned"] = self.get_nx_binned()
        O["ny_binned"] = self.get_ny_binned()
        O["pixel_size_binned"] = self.get_pixel_size("binned")
        O["pixel_size_unbinned"] = self.get_pixel_size("unbinned")
        O["distance"] = self.distance
        return O

    def _get_c(self,c_in,option):
        # Take into account binning if required
        if option == 'unbinned':
            pass
        elif option == 'binned':
            if (c_in % 1) == 0.5:
                c_out = (c_in-(self._binning-1)/2.)/(1.0*self._binning)
            else:
                c_out = c_in/(1.0*self._binning)
        else:
            logger.error("No valid option chosen.")
        return c_out

    def get_minimum_center_edge_distance(self,cx,cy):
        icx = self._nx-cx
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
    
    def get_absqx_max(self,wavelength,cx_unbinned):
        x_max = max([self._get_c(cx_unbinned,'binned'),self._mask.shape[1]-1-self._get_c(cx_unbinned,'binned')]) * self.get_pixel_size('binned')
        R_Ewald = 2*numpy.pi/wavelength
        phi = numpy.arctan2(x_max,self.distance)
        return 2 * R_Ewald * numpy.sin(phi/2.0)

    def get_absqy_max(self,wavelength,cy_unbinned):
        y_max = max([self._get_c(cy_unbinned,'binned'),self._mask.shape[0]-1-self._get_c(cy_unbinned,'binned')]) * self.get_pixel_size('binned')
        R_Ewald = 2*numpy.pi/wavelength
        phi = numpy.arctan2(y_max,self.distance)
        return 2 * R_Ewald * numpy.sin(phi/2.0)

    def get_absqz_max(self,wavelength,cx_unbinned,cy_unbinned):
        absqx_max = self.get_absqx_max(wavelength,cx_unbinned)
        absqy_max = self.get_absqy_max(wavelength,cy_unbinned)
        R_Ewald = 2*numpy.pi/wavelength
        phi = numpy.arcsin(numpy.sqrt(absqx_max**2+absqy_max**2)/R_Ewald)
        return R_Ewald * (1-numpy.cos(phi))

    def get_absq_max(self,wavelength,cx_unbinned,cy_unbinned):
        return max([self.get_absqx_max(wavelength,cx_unbinned),self.get_absqy_max(wavelength,cy_unbinned)])

    def get_real_space_resolution_element(self,wavelength,cx_unbinned,cy_unbinned):
        dX = numpy.pi / self.get_absq_max(wavelength,cx_unbinned,cy_unbinned)
        return dX

    def get_real_space_resolution_element_min(self,wavelength,cx_unbinned,cy_unbinned):
        if self._center_spread_limit == 0:
            return self.get_real_space_resolution_element(wavelength,cx_unbinned,cy_unbinned)
        else:
            dx,dy = numpy.meshgrid(range(3),range(3))
            dx -= 1
            dy -= 1
            dX = (numpy.pi / self.get_absq_max(wavelength,cx_unbinned+self._center_spread_limit/2.*dx,cy_unbinned+self._center_spread_limit/2.*dy)).min()
        return dX

    def get_max_achievable_crystallographic_resolution(self,wavelength,cx_unbinned,cy_unbinned):
        dx = 2 * numpy.pi / self.get_absqx_max(wavelength,cx_unbinned)
        dy = 2 * numpy.pi / self.get_absqy_max(wavelength,cy_unbinned)
        return [dx,dy]

    def detect_photons(self,I):
        I_det = self._noise.get(I)
        if self.saturation_level is not None:
            temp = I_det > self.saturation_level
            if temp.sum() > 0:
                I_det[temp] = self.saturation_level
        return I_det

