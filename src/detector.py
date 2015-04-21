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

import sys,os
import h5py
sys.path.append("utils")
import numpy

import logging
logger = logging.getLogger("Condor")
import utils.log
from utils.log import log 

import config
from utils.variation import Variation
from utils.pixelmask import PixelMask


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
        - pixel_size: Edge length of square pixel.
        - cx: Horizontal beam position in pixel. If argument is \'None\' or not given center is set to the middle. [None]
        - cy: Vertical beam position in pixel. If argument is \'None\' or not given center is set to the middle. [None]
        - center_variation: Variation of the center position. Either None, \'normal\', \'uniform\' or \'range\'. [None]
        - center_spread_x: If \'normal\', \'uniform\' or \'range\' center variation is specified spread defines width of the distribution in x.
        - center_spread_y: If \'normal\', \'uniform\' or \'range\' center variation is specified spread defines width of the distribution in y.
        - center_variation_n: I \'range\' center variation is specified this argument specifies the number of samples within the specified range.
        - noise: Noise that is put on the intensities when read. Either None, \'poisson\', \'normal\', \'uniform\' or \'normal_poisson or \'file\' or \'file_poisson\''. [None]
        - noise_spread: If \'normal\', \'uniform\' or \'normal_poisson\' noise is specified spread defines width of gaussian/uniform distribution.
        - noise_filename: If \'file\' or \'file_poisson\' noise is specified this HDF5 file contains the dataset that is added to an image.
        - noise_dataset:  If \'file\' or \'file_poisson\' noise is specified this HDF5 file dataset contains the signal that is added to a diffraction image. If the dataset has 3 dimensions a random frame (i.e. dataset[i_random,:,:]) is picked for every read out.
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
          + True: (pixels & PixelMask.PIXEL_IS_IN_MASK_DEFAULT) == 0 are the valid pixels
        
          EITHER
          - mask: 2d mask array. In the array 1 stands for valid pixel and 0 for an invalid pixel.
          OR:
          - mask_filename: path to HDF5 file where the mask is stored
          - mask_dataset: dataset path in file
        """

        # Check for valid set of keyword arguments
        req_keys = ["distance","pixel_size","cx","cy"]
        opt_keys = ["x_gap_size_in_pixel","y_gap_size_in_pixel","hole_diameter_in_pixel",
                    "noise","noise_spread","noise_variation_n",
                    "noise_filename","noise_dataset",
                    "center_variation","center_spread_x","center_spread_y","center_variation_n","center_spread_limit",
                    "saturation_level",
                    "mask","mask_filename","mask_dataset",
                    "nx","ny","downsampling"]
        miss_keys,ill_keys = config.check_input(kwargs.keys(),req_keys,opt_keys)
        if len(miss_keys) > 0: 
            for k in miss_keys:
                log(logger.error,"Cannot initialize Detector instance. %s is a necessary keyword." % k)
            sys.exit(1)
        if len(ill_keys) > 0:
            for k in ill_keys:
                log(logger.error,"Cannot initialize Detector instance. %s is an illegal keyword." % k)
            sys.exit(1)

        # Start initialisation            
        self.distance = kwargs["distance"]
        self.pixel_size = kwargs["pixel_size"]
        gx = kwargs["x_gap_size_in_pixel"]
        gy = kwargs["y_gap_size_in_pixel"]
        hd = kwargs["hole_diameter_in_pixel"]
        self.set_noise(noise=kwargs["noise"],
                       noise_spread=kwargs.get("noise_spread",None),
                       noise_variation_n=kwargs.get("noise_variation_n",None),
                       noise_filename=kwargs.get("noise_filename",None),
                       noise_dataset=kwargs.get("noise_dataset",None))
        self.set_center_variation(center_variation=kwargs["center_variation"],
                                  center_spread_x=kwargs.get("center_spread_x",None),
                                  center_spread_y=kwargs.get("center_spread_y",None),
                                  center_variation_n=kwargs.get("center_variation_n",None),
                                  center_spread_limit=kwargs.get("center_spread_limit",None))
        self.saturation_level = kwargs["saturation_level"]
        if "mask" in kwargs or  ("mask_filename" in kwargs and "mask_dataset" in kwargs):
            if "mask" in kwargs:
                M = kwargs["mask"]
            else:
                with h5py.File(kwargs["mask_filename"],"r") as f:
                    M = f[kwargs["mask_dataset"]][:,:]
            if not kwargs.get("mask_CXI_bitmask",True): M = (M & PixelMask.PIXEL_IS_IN_MASK_DEFAULT) == 0
            self.init_mask(mask=M)
        else:
            reqk = ["nx","ny"]
            for k in reqk:
                if k not in kwargs.keys():
                    log(logger.error,"Cannot initialize Detector instance. %s is a necessary keyword if no mask is given." % k)
                    sys.exit(1)
            self.init_mask(nx=kwargs["nx"],ny=kwargs["ny"],
                           x_gap_size_in_pixel=gx,y_gap_size_in_pixel=gy,hole_diameter_in_pixel=hd,
                           cx_hole=kwargs.get("cx_hole","middle"),cy_hole=kwargs.get("cy_hole","middle"))
        self.cx_mean = kwargs.get("cx","middle")
        self.cy_mean = kwargs.get("cy","middle")
        self.downsampling = kwargs.get("downsampling",None)

    def set_noise(self, noise, noise_spread, noise_variation_n, noise_filename, noise_dataset):
        if noise in ["file","file_poisson"]:
            self._noise_filename = noise_filename
            self._noise_dataset = noise_dataset
            self._noise = Variation("poisson" if noise == "file_poisson" else None,noise_spread,noise_variation_n,number_of_dimensions=1,name="noise")
        else:
            self._noise_filename = None
            self._noise_dataset = None
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
        - mask: array that defines the mask (CXI bitmask, uint16)
        
        OR:
        - nx: horizontal dimension of the mask (inside the mask, without counting any gaps)
        - ny: vertical dimension of the mask (inside the mask, without counting any gaps)
        - x_gap_size_in_pixel: horizontal gap is generated with given width (in unit pixel)
        - y_gap_size_in_pixel: vertical gap is generated with given width (in unit pixel)
        - hole_diameter_in_pixel: holeis generated with given diameter (in unit pixel)
        - cx_hole: horizontal center coordinate of the hole
        - cy_hole: vertical center coordinate of the hole

        """

        # init mask array
        if kwargs.get('mask',None) != None: 
            self._mask = numpy.array(kwargs['mask'],dtype=numpy.uint16)
            self._nx = self._mask.shape[1]
            self._ny = self._mask.shape[0]
        elif 'nx' in kwargs and 'ny' in kwargs:
            self._nx = kwargs["nx"]
            self._ny = kwargs["ny"]
            Nx = kwargs['nx']
            Ny = kwargs['ny']
            if 'x_gap_size_in_pixel' in kwargs:
                x_gap = kwargs['x_gap_size_in_pixel']
                Ny += x_gap
            if 'y_gap_size_in_pixel' in kwargs: 
                y_gap = kwargs['y_gap_size_in_pixel']
                Nx += y_gap
            self._mask = numpy.zeros(shape=(Ny,Nx),dtype=numpy.uint16)
        else:
            log(logger.error,"Either \'mask_array\' or \'nx\' and \'ny\' have to be specified.")
            return 
        self._nx = self._mask.shape[1]
        self._ny = self._mask.shape[0]

        # set pixels in gap to zero
        if 'x_gap_size_in_pixel' in kwargs:
            if kwargs['x_gap_size_in_pixel'] != 0:
                cy = numpy.ceil((self._ny-1)/2.)
                self._mask[cy-kwargs['x_gap_size_in_pixel']/2:cy-kwargs['x_gap_size_in_pixel']/2+kwargs['x_gap_size_in_pixel'],:] |= PixelMask.PIXEL_IS_MISSING
        if 'y_gap_size_in_pixel' in kwargs:
            if kwargs['y_gap_size_in_pixel'] != 0:
                cx = numpy.ceil((self._nx-1)/2.)
                self._mask[:,cx-kwargs['y_gap_size_in_pixel']/2:cx-kwargs['y_gap_size_in_pixel']/2+kwargs['y_gap_size_in_pixel']] |= PixelMask.PIXEL_IS_MISSING
        if 'hole_diameter_in_pixel' in kwargs:
            if kwargs['hole_diameter_in_pixel'] != 0:
                cx_hole = kwargs.get("cx_hole","middle")
                if cx_hole == "middle":
                    cx_hole = (self._nx-1)/2.
                cy_hole = kwargs.get("cy_hole","middle")
                if cy_hole == "middle":
                    cy_hole = (self._ny-1)/2.
                X,Y = numpy.meshgrid(numpy.arange(0,self._nx,1.0),
                                     numpy.arange(0,self._ny,1.0))
                X = X-cx_hole
                Y = Y-cy_hole
                R = numpy.sqrt(X**2 + Y**2)
                self._mask[R<=kwargs['hole_diameter_in_pixel']/2.0] |= PixelMask.PIXEL_IS_MISSING
        

    def get_mask(self,intensities,output_bitmask=True):
        if output_bitmask:
            return self.get_bitmask(intensities)
        else:
            return numpy.array(self.get_bitmask(intensities) == 0,dtype="bool")
    
    def get_bitmask(self,intensities):
        M = self._mask.copy()
        if self.saturation_level is not None:
            M[intensities >= self.saturation_level] |= PixelMask.PIXEL_IS_SATURATED
        return M

    def get_cx_mean_value(self):
        if self.cx_mean == "middle":
            return (self._nx-1) / 2.
        else:
            return self.cx_mean

    def get_cy_mean_value(self):
        if self.cy_mean == "middle":
            return (self._ny-1) / 2.
        else:
            return self.cy_mean

    def get_next(self):
        O = {}
        cx_mean = self.get_cx_mean_value()
        cy_mean = self.get_cy_mean_value()
        # Center variation
        ready = False 
        while not ready:
            cy_now, cx_now = self._center_variation.get([cy_mean,cx_mean])
            ready = True
            if self._center_spread_limit > 0:
                if (self._center_spread_limit < abs(cx_now - cx_mean)*2) and (self._center_spread_limit < abs(cy_now - cy_mean)*2):
                    ready = False
        O["cx"] = cx_now
        O["cy"] = cy_now
        O["solid_angle_pixel"] = self.get_pixel_solid_angle()
        O["nx"] = self._nx
        O["ny"] = self._ny
        O["pixel_size"] = self.pixel_size
        O["distance"] = self.distance
        if self.downsampling is not None:
            O["cx_xxx"] = utils.resample.downsample_pos(cx_now,self._nx,self.downsampling)
            O["cy_xxx"] = utils.resample.downsample_pos(cy_now,self._ny,self.downsampling)
        return O

    def get_minimum_center_edge_distance(self,cx,cy):
        icx = self._nx-cx
        icy = self._ny-cy
        return min([cx,icx,cy,icy])*self.pixel_size

    def get_pixel_solid_angle(self):
        return self.pixel_size**2 / self.distance**2
    
    def get_absqx_max(self,wavelength,cx):
        x_max = max([cx,self._nx-1-cx]) * self.pixel_size
        R_Ewald = 2*numpy.pi/wavelength
        phi = numpy.arctan2(x_max,self.distance)
        return 2 * R_Ewald * numpy.sin(phi/2.0)

    def get_absqy_max(self,wavelength,cy):
        y_max = max([cy,self._ny-1-cy]) * self.pixel_size
        R_Ewald = 2*numpy.pi/wavelength
        phi = numpy.arctan2(y_max,self.distance)
        return 2 * R_Ewald * numpy.sin(phi/2.0)

    def get_absqz_max(self,wavelength,cx,cy):
        absqx_max = self.get_absqx_max(wavelength,cx)
        absqy_max = self.get_absqy_max(wavelength,cy)
        R_Ewald = 2*numpy.pi/wavelength
        phi = numpy.arcsin(numpy.sqrt(absqx_max**2+absqy_max**2)/R_Ewald)
        return R_Ewald * (1-numpy.cos(phi))

    def get_absq_max(self,wavelength,cx,cy):
        return max([self.get_absqx_max(wavelength,cx),self.get_absqy_max(wavelength,cy)])

    def get_real_space_resolution_element(self,wavelength,cx,cy):
        dX = numpy.pi / self.get_absq_max(wavelength,cx,cy)
        return dX

    def get_real_space_resolution_element_min(self,wavelength,cx,cy):
        if self._center_spread_limit == 0:
            return self.get_real_space_resolution_element(wavelength,cx,cy)
        else:
            dx,dy = numpy.meshgrid(range(3),range(3))
            dx -= 1
            dy -= 1
            cx_min = cx - self._center_spread_limit
            cx_max = cx + self._center_spread_limit
            cy_min = cy - self._center_spread_limit
            cy_max = cy + self._center_spread_limit
            cx_off = [cx_min,cx_max][numpy.array([self._nx-cx_min,cx_max]).argmax()]
            cy_off = [cy_min,cy_max][numpy.array([self._ny-cy_min,cy_max]).argmax()]
            return self.get_real_space_resolution_element(wavelength,cx_off,cy_off)

    def get_max_achievable_crystallographic_resolution(self,wavelength,cx,cy):
        dx = 2 * numpy.pi / self.get_absqx_max(wavelength,cx)
        dy = 2 * numpy.pi / self.get_absqy_max(wavelength,cy)
        return [dx,dy]

    def detect_photons(self,I):
        I_det = self._noise.get(I)
        if self._noise_filename is not None:
            with h5py.File(self._noise_filename,"r") as f:
                ds = f[self._noise_dataset]
                if len(list(ds.shape)) == 2:
                    bg = ds[:,:]
                else:
                    bg = ds[numpy.random.randint(ds.shape[0]),:,:]
            I_det = I_det + bg
        if self.saturation_level is not None:
            I_det = numpy.clip(I_det, -numpy.inf, self.saturation_level)
        M_det = self.get_bitmask(I_det)
        if self.downsampling is not None:
            IXxX_det, MXxX_det = utils.resample.downsample(I_det,self.downsampling,mode="integrate",
                                                           mask2d0=M_det,bad_bits=PixelMask.PIXEL_IS_IN_MASK,min_N_pixels=1)
        else:
            IXxX_det = None
            MXxX_det = None
        return I_det, M_det, IXxX_det, MXxX_det
