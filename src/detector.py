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
logger = logging.getLogger(__name__)

import condor.utils.log
from condor.utils.log import log_and_raise_error,log_warning,log_info,log_debug
from condor.utils.config import load_config
import utils.resample
from condor.utils.variation import Variation
from condor.utils.pixelmask import PixelMask
from condor.utils.linalg import length


def load_detector(conf=None):
    """
    Create new Detector instance and load parameters from a Condor configuration file or dictionary.
    
    Kwargs:
       :conf(str): Condor configuration file or dictionary (default = None)
    """
    C = condor.utils.config.load_config({"detector": load_config(conf)["detector"]}, {"detector": load_config(condor.CONDOR_default_conf)["detector"]})
    detector = Detector(**C["detector"])
    return detector

class Detector:
    """
    Class for photon area-detector
    """
    def __init__(self, distance, pixel_size,
                 x_gap_size_in_pixel=0, y_gap_size_in_pixel=0, hole_diameter_in_pixel=0, cx_hole="middle", cy_hole="middle",
                 noise=None, noise_spread=None, noise_variation_n=None, noise_filename=None, noise_dataset=None,
                 cx="middle", cy="middle", center_variation=None, center_spread_x=None, center_spread_y=None, center_variation_n=None, center_spread_limit=None,
                 saturation_level=None, mask=None, mask_filename=None, mask_dataset=None, mask_is_cxi_bitmask=False,
                 nx=None, ny=None, binning=None):
        """
        Initialisation of a Detector instance
        
        Args:
           :distance(float): Distance between detector and interaction point
           :pixel_size(float): Edge length of square pixel

        Kwargs:
           :cx: Horizontal beam position in pixel (fast changing dimension). If argument is None or not given center is set to the middle (default = None)
           :cy: Vertical beam position in pixel (slowly changing dimension). If argument is None or not given center is set to the middle (default = None)
           :center_variation(str): Variation of the center position. Either None, \"normal\", \"uniform\" or \"range\" (default = None)
           :center_spread_x(float): If \"normal\", \"uniform\" or \"range\" center variation is specified spread defines width of the distribution in x
           :center_spread_y(float): If \"normal\", \"uniform\" or \"range\" center variation is specified spread defines width of the distribution in y
           :center_variation_n(int): I \"range\" center variation is specified this argument specifies the number of samples within the specified range
           :noise(str): Noise that is put on the intensities when read. Either None, \"poisson\", \"normal\", \"uniform\" or \"normal_poisson or \"file\" or \"file_poisson\" (default = None)
           :noise_spread(float): If \"normal\", \"uniform\" or \"normal_poisson\" noise is specified spread defines width of gaussian/uniform distribution (default = None)
           :noise_filename(str): If \"file\" or \"file_poisson\" noise is specified this HDF5 file contains the dataset that is added to an image (default = None)
           :noise_dataset(str):  If \"file\" or \"file_poisson\" noise is specified this HDF5 file dataset contains the signal that is added to a diffraction image. If the dataset has 3 dimensions a random frame (i.e. dataset[i_random,:,:]) is picked for every read out (default = None)
           :saturation_level(float): Value at which detector pixels satutrate (default = None)
           :binning(int): Pixel binning factor, intensies are integrated over patches of \'binning\' x \'binning\' pixels (default = None)

        The mask can be specified in three alternative ways:

        By parametrisation
           :nx(int): Number of pixels in horizontal direction not including a potential gap (default = None)
           :ny(int): Number of pixels in vertical direction not including a potential gap (default = None)
           :x_gap_size_in_pixel(int): Width of central horizontal gap in pixel (default = 0)
           :y_gap_size_in_pixel(int): Width of central vertical gap in pixel (default = 0)
           :hole_diameter_in_pixel(int): Diameter of central hole in pixel (default = 0)
        
        or by specifying the mask with an array
           :mask_filename: path to HDF5 file where the mask is stored (default = None)
           :mask_dataset: dataset path in file (default = None)
           :mask_CXI_bitmask(bool): Bool wheter or not data is stored as a CXI bitmask. False: pixels with value 0 are invalid, pixels with value 1 are valid. True: (pixels & PixelMask.PIXEL_IS_IN_MASK_DEFAULT) == 0 are the valid pixels (default = False)

        or by reading the mask from a dataset in an HDF5 file
           :mask: 2d mask as numpy array (default = None)
           :mask_CXI_bitmask(bool): Bool wheter or not data is stored as a CXI bitmask. False: pixels with value 0 are invalid, pixels with value 1 are valid. True: (pixels & PixelMask.PIXEL_IS_IN_MASK_DEFAULT) == 0 are the valid pixels (default = False)
        """

        self.distance = distance
        self.pixel_size = pixel_size
        self._init_mask(mask=mask, mask_is_cxi_bitmask=mask_is_cxi_bitmask, mask_filename=mask_filename, mask_dataset=mask_dataset, nx=nx, ny=ny,
                        x_gap_size_in_pixel=x_gap_size_in_pixel, y_gap_size_in_pixel=y_gap_size_in_pixel, cx_hole=cx_hole, cy_hole=cy_hole, hole_diameter_in_pixel=hole_diameter_in_pixel)
        self.cx_mean = cx
        self.cy_mean = cy
        self.set_center_variation(center_variation=center_variation,
                                  center_spread_x=center_spread_x,
                                  center_spread_y=center_spread_y,
                                  center_variation_n=center_variation_n,
                                  center_spread_limit=center_spread_limit)
        self.set_noise(noise=noise,
                       noise_spread=noise_spread,
                       noise_variation_n=noise_variation_n,
                       noise_filename=noise_filename,
                       noise_dataset=noise_dataset)
        self.saturation_level = saturation_level
        self.binning = binning

    def get_conf(self):
        conf = {}
        conf["detector"] = {}
        conf["detector"]["distance"]           = self.distance
        conf["detector"]["pixel_size"]         = self.pixel_size
        conf["detector"]["cx"]                 = self.cx_mean
        conf["detector"]["cy"]                 = self.cy_mean
        cvar = self._center_variation.get_conf()
        conf["detector"]["center_variation"]   = cvar["mode"]
        conf["detector"]["center_spread_x"]    = cvar["spread"][1]
        conf["detector"]["center_spread_y"]    = cvar["spread"][0]
        conf["detector"]["center_variation_n"] = cvar["n"]
        noise = self._noise.get_conf()
        conf["detector"]["noise"]              = noise["mode"]
        conf["detector"]["noise_spread"]       = noise["spread"]
        conf["detector"]["noise_filename"]     = self._noise_filename
        conf["detector"]["noise_dataset"]      = self._noise_dataset
        conf["detector"]["saturation_level"]   = self.saturation_level
        conf["detector"]["mask"]               = self._mask.copy()
        conf["detector"]["mask_CXI_bitmask"]   = True
        return conf
        
    def set_noise(self, noise=None, noise_spread=None, noise_variation_n=None, noise_filename=None, noise_dataset=None):
        """
        Set detector noise type and parameters

        Kwargs:
           :noise(str): Noise that is put on the intensities when read. Either None, \"poisson\", \"normal\", \"uniform\" or \"normal_poisson or \"file\" or \"file_poisson\" (default = None)
           :noise_spread(float): If \"normal\", \"uniform\" or \"normal_poisson\" noise is specified spread defines width of gaussian/uniform distribution (default = None)
           :noise_filename(str): Location of HDF5 file that contains noise data that is added to simulated data. Specify the dataset with the argument \'noise_dataset\' (see below). The arguments \'noise_filename\' and \'noise_dataset\' only matter if noise=\"file\" or noise=\"file_poisson\" (default = None)
           :noise_dataset(str): Dataset that contains the noise data that is added to simulated data. If the dataset has 3 dimensions a random frame (i.e. dataset[i_random,:,:]) is picked for every simulated pattern. Specify the filename with the argument \'noise_filename\' (see above). The arguments \'noise_filename\' and \'noise_dataset\' only matter if noise=\"file\" or noise=\"file_poisson\"  (default = None)
        """
        if noise in ["file","file_poisson"]:
            self._noise_filename = noise_filename
            self._noise_dataset = noise_dataset
            self._noise = Variation("poisson" if noise == "file_poisson" else None, noise_spread, noise_variation_n, number_of_dimensions=1, name="noise")
        else:
            self._noise_filename = None
            self._noise_dataset = None
            self._noise = Variation(noise, noise_spread, noise_variation_n, number_of_dimensions=1, name="noise")

    def set_center_variation(self, center_variation, center_spread_x, center_spread_y, center_variation_n, center_spread_limit):
        self._center_variation = Variation(center_variation,[center_spread_y,center_spread_x],center_variation_n,number_of_dimensions=2,name="center position")
        self._center_spread_limit = center_spread_limit
        
    def _init_mask(self, mask, mask_is_cxi_bitmask, mask_filename, mask_dataset, nx, ny, x_gap_size_in_pixel, y_gap_size_in_pixel, cx_hole, cy_hole, hole_diameter_in_pixel):
        if mask is not None or (mask_filename is not None and mask_dataset is not None):
            if mask is not None:
                # Copy mask from array
                self._mask = numpy.array(mask, dtype=numpy.uint16)
            else:
                # Read mask from file
                with h5py.File(mask_filename,"r") as f:
                    self._mask = numpy.array(f[mask_dataset][:,:], dtype=numpy.uint16)
            if not mask_is_cxi_bitmask:
                # Convert maskt to CXI bit format
                self._mask = (self._mask == 0) * PixelMask.PIXEL_IS_MISSING
        elif nx is not None and ny is not None:
            # Initialise empty mask
            self._mask = numpy.zeros(shape=(ny+y_gap_size_in_pixel, nx+x_gap_size_in_pixel),dtype=numpy.uint16)
        else:
            log_and_raise_error(logger, "Either \"mask\" or \"nx\" and \"ny\" have to be specified.")
            sys.exit(1)
        self._nx = self._mask.shape[1]
        self._ny = self._mask.shape[0]
        # Mask out pixels in gaps
        if x_gap_size_in_pixel > 0:
            cy = numpy.ceil((self._ny-1)/2.)
            self._mask[cy-x_gap_size_in_pixel/2:cy-x_gap_size_in_pixel/2+x_gap_size_in_pixel,:] |= PixelMask.PIXEL_IS_MISSING
        if y_gap_size_in_pixel > 0:
            cx = numpy.ceil((self._nx-1)/2.)
            self._mask[:,cx-y_gap_size_in_pixel/2:cx-y_gap_size_in_pixel/2+y_gap_size_in_pixel] |= PixelMask.PIXEL_IS_MISSING
        # Mask out pixels in hole    
        if hole_diameter_in_pixel > 0:
            if cx_hole == "middle":
                cx_hole = (self._nx-1)/2.
            if cy_hole == "middle":
                cy_hole = (self._ny-1)/2.
            Y,X = numpy.indices((self._ny,self._nx), dtype=numpy.float64)
            X = X-cx_hole
            Y = Y-cy_hole
            R = numpy.sqrt(X**2 + Y**2)
            tmp = R<=hole_diameter_in_pixel/2.0
            if tmp.sum() > 0:
                self._mask[tmp] |= PixelMask.PIXEL_IS_MISSING
        
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

    def _get_c_mean_value(self):
        c =  numpy.array([self.get_cy_mean_value(), self.get_cx_mean_value()])
        return c
    
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
            if self._center_variation._mode is not None:
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
        if self.binning is not None:
            O["cx_xxx"] = utils.resample.downsample_pos(cx_now,self._nx,self.binning)
            O["cy_xxx"] = utils.resample.downsample_pos(cy_now,self._ny,self.binning)
        return O

    def get_pixel_solid_angle(self):
        omega = self.pixel_size**2 / self.distance**2
        return omega

    def get_x_max(self, cx = None, cy = None, center_variation = False):
        if cx is None or center_variation:
            cx = self.get_cx_mean_value()
        lim = self._center_spread_limit if (center_variation and self._center_spread_limit is not None) else 0
        cx_min = cx - lim/2.
        cx_max = cx + lim/2.
        dist_max = max([self._nx-1-cx_min, cx_max]) * self.pixel_size
        return dist_max

    def get_y_max(self, cx = None, cy = None, center_variation = False):
        if cy is None or center_variation:
            cy = self.get_cy_mean_value()
        lim = self._center_spread_limit if (center_variation and self._center_spread_limit is not None) else 0
        cy_min = cy - lim/2.
        cy_max = cy + lim/2.
        dist_max = max([self._ny-1-cy_min, cy_max]) * self.pixel_size
        return dist_max

    def get_p_max(self, cx = None, cy = None, pos = "corner", center_variation = False):
        x_max = self.get_x_max(cx=cx, cy=cy, center_variation=center_variation)
        y_max = self.get_y_max(cx=cx, cy=cy, center_variation=center_variation)
        log_debug(logger, "y_max = %.1f pix, x_max = %.1f pix" % (y_max/self.pixel_size,x_max/self.pixel_size))
        if pos == "corner":
            return numpy.array([self.distance, y_max, x_max])
        elif pos == "edge":
            if x_max > y_max:
                return numpy.array([self.distance, 0., x_max])
            else:
                return numpy.array([self.distance, y_max, 0.])
        else:
            log_and_raise_error(logger, "Invalid input: pos=%s. Input must be either \"corner\" or \"edge\"." % pos)

    def get_q_max(self, wavelength, cx = None, cy = None, pos = "corner", center_variation = False):
        p = self.get_p_max(cx=cx, cy=cy, pos=pos, center_variation=center_variation)
        q = q_from_p(p, wavelength)
        return q

    def _get_resolution_element(self, wavelength, cx = None, cy = None, center_variation = False, pos="edge"):
        res = numpy.pi / self.get_q_max(wavelength, cx=cx, cy=cy, pos=pos, center_variation=center_variation)
        return res

    def get_resolution_element_z(self, wavelength, cx = None, cy = None, center_variation = False):
        return self._get_resolution_element(wavelength, cx=cx, cy=cy, center_variation=center_variation, pos="corner")[0]

    def get_resolution_element_y(self, wavelength, cx = None, cy = None, center_variation = False):
        return self._get_resolution_element(wavelength, cx=cx, cy=cy, center_variation=center_variation, pos="corner")[1]

    def get_resolution_element_x(self, wavelength, cx = None, cy = None, center_variation = False):
        return self._get_resolution_element(wavelength, cx=cx, cy=cy, center_variation=center_variation, pos="corner")[2]

    def get_resolution_element_r(self, wavelength, cx = None, cy = None, center_variation = False):
        qmax = self.get_q_max(wavelength, cx=cx, cy=cy, pos="corner", center_variation=center_variation)
        res = numpy.pi / length(qmax)
        return res
        
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
        if self.binning is not None:
            IXxX_det, MXxX_det = utils.resample.downsample(I_det,self.binning,mode="integrate",
                                                           mask2d0=M_det,bad_bits=PixelMask.PIXEL_IS_IN_MASK,min_N_pixels=1)
        else:
            IXxX_det = None
            MXxX_det = None
        return I_det, M_det, IXxX_det, MXxX_det

def q_from_p(p, wavelength):
    p0 = p / length(p)
    R_Ewald = 2*numpy.pi / wavelength
    k0 = R_Ewald * numpy.array([1.,0.,0.])
    k1 = R_Ewald * p0
    q = k0 - k1
    return q
