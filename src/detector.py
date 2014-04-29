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
        - binning: Number of binned pixels (binning x binning). [1]
        - noise: Noise that is put on the intensities when read. Either \'none\' or \'poisson\'. [\'none\']
        - parent: Input object that includes detector object. This variable is optional. [None]

        EITHER (default):
        - nx: Number of pixels in horizontal direction not including a potential gap.
        - Ny: Number of pixels in vertical direction not including a potential gap.
        - x_gap_size_in_pixel: Width of central horizontal gap in pixel. [0]
        - y_gap_size_in_pixel: Width of central vertical gap in pixel. [0]
        - hole_diameter_in_pixel: Diameter of central hole in pixel. [0]
        
        OR:
        - mask_CXI_bitmask: bool wheter or not data is stored as a CXI bitmask. [False]
          + False: pixels with value 0 are invalid, pixels with value 1 are valid.
          + True: (pixels & cxitools.PIXEL_IS_IN_MASK) == 0 are the valid pixels
        
          EITHER
          - mask: 2d mask array. In the array 1 stands for valid pixel and 0 for an invalid pixel.
          OR:
          - mask_filename: path to HDF5 file where the mask is stored
          - mask_dataset: dataset path in file
        """

        reqk = ["distance","pixel_size"]
        for k in reqk:
            if k not in kwargs.keys():
                logger.error("Cannot initialize Detector instance. %s is a necessary keyword." % k)
                return

        self._parent = kwargs.get("parent",None)
        self.distance = kwargs["distance"]
        self.pixel_size = kwargs["pixel_size"]
        gx = kwargs.get('x_gap_size_in_pixel',0)
        gy = kwargs.get('y_gap_size_in_pixel',0)
        hd = kwargs.get('hole_diameter_in_pixel',0)
        self.noise = kwargs.get('noise','none')
        self.binning = kwargs.get('binning',1)
        if 'mask' in kwargs or ('mask_filename' in kwargs.keys() and "mask_dataset" in kwargs.keys()):
            if "mask" in kwargs:
                M = kwargs["mask"]
            else:
                f = h5py.File(kwargs["mask_filename"],"r")
                M = f["mask_dataset"][:]
                f.close()
            if kwargs.get("mask_CXI_bitmask",False): M = (M & cxitools.PIXEL_IS_IN_MASK) == 0
            self.init_mask(mask=M)
            self.cx = kwargs.get('cx',(self.mask.shape[1]-1)/(2.*self.binning))
            self.cy = kwargs.get('cy',(self.mask.shape[0]-1)/(2.*self.binning))              
        else:
            reqk = ["nx","ny"]
            for k in reqk:
                if k not in kwargs.keys():
                    logger.error("Cannot initialize Detector instance. %s is a necessary keyword if no mask is given." % k)
                    return
            self.Nx = kwargs["nx"]
            self.Ny = kwargs["ny"]
            self.cx = kwargs.get('cx',(self.Nx+gy-1)/2.)
            self.cy = kwargs.get('cy',(self.Ny+gx-1)/2.)   
            self.init_mask(nx=kwargs["nx"],ny=kwargs["ny"],x_gap_size_in_pixel=gx,y_gap_size_in_pixel=gy,hole_diameter_in_pixel=hd,binning=self.binning)

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
            self.binning = kwargs['binning']
        else:
            self.binning = 1
        if kwargs.get('mask',None) != None: 
            self.mask = kwargs['mask']
            self.Nx = self.mask.shape[1]*self.binning
            self.Ny = self.mask.shape[0]*self.binning
        elif 'nx' in kwargs and 'ny' in kwargs:
            self.Nx = kwargs["nx"]
            self.Ny = kwargs["ny"]
            Nx = kwargs['nx']/self.binning; Ny = kwargs['ny']/self.binning
            if 'x_gap_size_in_pixel' in kwargs:
                x_gap = kwargs['x_gap_size_in_pixel']/self.binning
                Ny += x_gap
            if 'y_gap_size_in_pixel' in kwargs: 
                y_gap = kwargs['y_gap_size_in_pixel']/self.binning
                Nx += y_gap
            self.mask = numpy.ones(shape=(Ny,Nx))
        else:
            logger.error("Either \'mask_array\' or \'nx\' and \'ny\' have to be specified.")
            return 
            
        # set pixels in gap to zero
        if 'x_gap_size_in_pixel' in kwargs:
            if kwargs['x_gap_size_in_pixel'] != 0:
                cy = numpy.ceil((self.mask.shape[0]-1)/2.)
                self.mask[cy-kwargs['x_gap_size_in_pixel']/2/self.binning:cy-kwargs['x_gap_size_in_pixel']/2/self.binning+kwargs['x_gap_size_in_pixel']/self.binning,:] = 0
        if 'y_gap_size_in_pixel' in kwargs:
            if kwargs['y_gap_size_in_pixel'] != 0:
                cx = numpy.ceil((self.mask.shape[1]-1)/2.)
                self.mask[:,cx-kwargs['y_gap_size_in_pixel']/2/self.binning:cx-kwargs['y_gap_size_in_pixel']/2/self.binning+kwargs['y_gap_size_in_pixel']/self.binning] = 0
        if 'hole_diameter_in_pixel' in kwargs:
            if kwargs['hole_diameter_in_pixel'] != 0:
                cy = (self.mask.shape[0]-1)/2.
                cx = (self.mask.shape[1]-1)/2.
                X,Y = numpy.meshgrid(numpy.arange(0,self.mask.shape[1],1.0),
                                     numpy.arange(0,self.mask.shape[0],1.0))
                X = X-cx
                Y = Y-cy
                R = numpy.sqrt(X**2 + Y**2)
                self.mask[R<=kwargs['hole_diameter_in_pixel']/(2.0*self.binning)] = 0

    def set_cy(self,cy=None):
        """
        Function sets vertical center position:
        =======================================
        Arguments:
        cy : vertical center position in pixel. If argument is None or not given center is set to the middle.
        """
        if cy == None:
            self.cy = (self.mask.shape[0]-1)/2.
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
            self.cx = (self.mask.shape[1]-1)/2.
        else:
            self.cx = cx

    def get_cx(self,option='unbinned'):
        if self.cx == "middle":
            cx = (self.Nx-1)/2.
        else:
            cx = self.cx
        if option == 'unbinned':
            return cx
        elif option == 'binned':
            if (cx % 1) == 0.5:
                return (cx-(self.binning-1)/2.)/(1.0*self.binning)
            else:
                return cx/(1.0*self.binning)
        else:
            logger.error("No valid option chosen.")

    def get_cy(self,option='unbinned'):
        if self.cy == "middle":
            cy = (self.Ny-1)/2.
        else:
            cy = self.cy
        if option == 'unbinned':
            return cy
        elif option == 'binned':
            if (cy % 1) == 0.5:
                return (cy-(self.binning-1)/2.)/(1.0*self.binning)
            else:
                return cy/(1.0*self.binning)
        else:
            logger.error("No valid option chosen.")
    def get_minimum_center_edge_distance(self):
        cx = self.get_cx()
        icx = self.Nx-cx
        cy = self.get_cy()
        icy = self.Ny-cy
        return min([cx,icx,cy,icy])*self.get_pixel_size()
    def get_pixel_size(self,option='unbinned'):
        if option == 'unbinned':
            return self.pixel_size
        elif option == 'binned':
            return self.pixel_size*self.binning
        else:
            logger.error("No valid option chosen.")

    def get_pixel_solid_angle(self,option='unbinned'):
        return self.get_pixel_size(option)**2 / self.distance**2
    
    def generate_absqmap(self,**kwargs):
        X,Y = numpy.meshgrid(numpy.arange(self.mask.shape[1]),
                             numpy.arange(self.mask.shape[0]))
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
        X,Y = numpy.meshgrid(numpy.arange(self.mask.shape[1]),
                             numpy.arange(self.mask.shape[0]))
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
            #qmap /= self.get_absq_max()/0.5
            qmap /= self.get_absq_max()/0.5*numpy.sqrt(2)
        return qmap

    def get_absqx_max(self,**kwargs):
        if "wavelength" in kwargs:
            w = kwargs["wavelength"]
        else:
            w = self._parent.source.photon.get_wavelength()
        x_max = max([self.get_cx('binned'),self.mask.shape[1]-1-self.get_cx('binned')]) * self.get_pixel_size('binned')
        R_Ewald = 2*numpy.pi/w
        phi = numpy.arctan2(x_max,self.distance)
        return 2 * R_Ewald * numpy.sin(phi/2.0)

    def get_absqy_max(self,**kwargs):
        if "wavelength" in kwargs:
            w = kwargs["wavelength"]
        else:
            w = self._parent.source.photon.get_wavelength()
        y_max = max([self.get_cy('binned'),self.mask.shape[0]-1-self.get_cy('binned')]) * self.get_pixel_size('binned')
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

    def detect_photons(self,data):
        if self.noise == "poisson":
            return numpy.random.poisson(data)
        else:
            return data
