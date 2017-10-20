# -----------------------------------------------------------------------------------------------------
# CONDOR
# Simulator for diffractive single-particle imaging experiments with X-ray lasers
# http://xfel.icm.uu.se/condor/
# -----------------------------------------------------------------------------------------------------
# Copyright 2016 Max Hantke, Filipe R.N.C. Maia, Tomas Ekeberg
# Condor is distributed under the terms of the BSD 2-Clause License
# -----------------------------------------------------------------------------------------------------
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# 
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# -----------------------------------------------------------------------------------------------------
# General note:
# All variables are in SI units by default. Exceptions explicit by variable name.
# -----------------------------------------------------------------------------------------------------

import sys,os
import collections
sys.path.append("utils")
import numpy

import logging
logger = logging.getLogger(__name__)

import condor.utils.log
from condor.utils.log import log_and_raise_error,log_warning,log_info,log_debug

import utils.resample
from condor.utils.variation import Variation
from condor.utils.pixelmask import PixelMask
from condor.utils.linalg import length
import condor.utils.testing
import condor.utils.scattering_vector


class Detector:
    """
    Class for a photon area-detector

    .. image:: images/detector_schematic.jpg
        
    **Arguments:**

      :distance (float): Distance from interaction point to detector plane

      :pixel_size (float): Edge length of detector pixel (square shape)

    **Keyword arguments:**
    
      :cx (float): Horizontal beam position in unit pixel. If ``cx=None`` beam will be positioned in the center (default ``None``)

      :cy (float): Vertical beam position in unit pixel If ``cy=None`` beam will be positioned in the center (default ``None``)

      :center_variation (str): See :meth:`condor.detector.Detector.set_center_variation` (default ``None``)
        
      :center_spread_x (float): See :meth:`condor.detector.Detector.set_center_variation` (default ``None``)
    
      :center_spread_y (float): See :meth:`condor.detector.Detector.set_center_variation` (default ``None``)

      :center_variation_n (int): See :meth:`condor.detector.Detector.set_center_variation` (default ``None``)

      :noise (str): See :meth:`condor.detector.Detector.set_noise` (default ``None``)

      :noise_spread (float): See :meth:`condor.detector.Detector.set_noise` (default ``None``)

      :noise_filename (str): See :meth:`condor.detector.Detector.set_noise` (default ``None``)

      :noise_dataset (str):  See :meth:`condor.detector.Detector.set_noise` (default ``None``)

      :saturation_level (float): Value at which detector pixels satutrate (default ``None``)

      :binning (int): Pixel binning factor, intensies are integrated over square patches that have an area of ``binning`` x ``binning`` pixels (default ``None``)

      :mask_CXI_bitmask (bool): If ``True`` the provided mask (``mask_dataset`` or ``mask``) is a CXI bitmask. For documentation on the implementation of CXI bitmasks see :class:`condor.utils.pixelmask.PixelMask` (default ``False``)

      :solid_angle_correction (bool): Whether or not solid angle correction shall be applied (default ``True``)

        *Choose one of the following options:*

        ==================== =============================================================================
        ``mask_CXI_bitmask`` valid pixels
        ==================== =============================================================================
        ``False``            ``1``
        ``True``             ``(pixels & condor.utils.pixelmask.PixelMask.PIXEL_IS_IN_MASK_DEFAULT) == 0``
        ==================== =============================================================================

      **There are 3 alternative options to specify shape and mask of the detector**

        *A) Parameters*

          :nx (int): Number of pixels in *x* direction (not including a potential gap or hole) (default ``None``)

          :ny (int): Number of pixels in *y* direction (not including a potential gap or hole) (default ``None``)

          :x_gap_size_in_pixel (int): Size of central gap along *x* in unit pixel (default ``None``)
    
          :y_gap_size_in_pixel (int): Size of central gap along *y* in unit pixel (default ``None``)

          :hole_diameter_in_pixel (int): Diameter of central hole in unit pixel (default ``None``)
    
        *B) HDF5 dataset for mask*

          :mask_filename (str): Location of HDF5 file that contains dataset for mask (default ``None``)

          :mask_dataset (str): HDF5 dataset (in the file specified by the argument ``mask_filename``) that contains the mask data. Toggle the option ``mask_CXI_bitmask`` for decoding options (default ``None``)


        *C) Numpy array for mask*

          :mask (array): 2D numpy integer array that defines the mask. Toggle ``mask_CXI_bitmask`` for decoding options (default ``None``)
    """
    def __init__(self, distance, pixel_size,
                 x_gap_size_in_pixel=0, y_gap_size_in_pixel=0, hole_diameter_in_pixel=0, cx_hole=None, cy_hole=None,
                 noise=None, noise_spread=None, noise_variation_n=None, noise_filename=None, noise_dataset=None,
                 cx=None, cy=None, center_variation=None, center_spread_x=None, center_spread_y=None, center_variation_n=None,
                 saturation_level=None, mask=None, mask_filename=None, mask_dataset=None, mask_is_cxi_bitmask=False, solid_angle_correction=True,
                 nx=None, ny=None, binning=None):

        self.distance = distance
        self.pixel_size = float(pixel_size)
        self._init_mask(mask=mask, mask_is_cxi_bitmask=mask_is_cxi_bitmask, mask_filename=mask_filename, mask_dataset=mask_dataset, nx=nx, ny=ny,
                        x_gap_size_in_pixel=x_gap_size_in_pixel, y_gap_size_in_pixel=y_gap_size_in_pixel, cx_hole=cx_hole, cy_hole=cy_hole, hole_diameter_in_pixel=hole_diameter_in_pixel)
        self.cx_mean = cx if cx != 'middle' else None
        self.cy_mean = cy if cy != 'middle' else None
        self.set_center_variation(center_variation=center_variation,
                                  center_spread_x=center_spread_x,
                                  center_spread_y=center_spread_y,
                                  center_variation_n=center_variation_n)
        self.set_noise(noise=noise,
                       noise_spread=noise_spread,
                       noise_variation_n=noise_variation_n,
                       noise_filename=noise_filename,
                       noise_dataset=noise_dataset)
        self.saturation_level = saturation_level
        self.binning = binning
        self.solid_angle_correction = solid_angle_correction

    def get_conf(self):
        """
        Get configuration in form of a dictionary. Another identically configured Detector instance can be initialised by:

        .. code-block:: python

          conf = D0.get_conf()         # D0: already existing Detector instance
          D1 = condor.Detector(**conf) # D1: new Detector instance with the same configuration as D0  
        """
        conf = {}
        conf["detector"] = {}
        conf["detector"]["distance"]           = self.distance
        conf["detector"]["pixel_size"]         = self.pixel_size
        conf["detector"]["cx"]                 = self.cx_mean
        conf["detector"]["cy"]                 = self.cy_mean
        cvar = self._center_variation.get_conf()
        conf["detector"]["center_variation"]   = cvar["mode"]
        conf["detector"]["center_spread_x"]    = cvar["spread"][0]
        conf["detector"]["center_spread_y"]    = cvar["spread"][1]
        conf["detector"]["center_variation_n"] = cvar["n"]
        noise = self._noise.get_conf()
        conf["detector"]["noise"]              = noise["mode"]
        conf["detector"]["noise_spread"]       = noise["spread"]
        conf["detector"]["noise_filename"]     = self._noise_filename
        conf["detector"]["noise_dataset"]      = self._noise_dataset
        conf["detector"]["saturation_level"]   = self.saturation_level
        conf["detector"]["mask"]               = self._mask.copy()
        conf["detector"]["mask_CXI_bitmask"]   = True
        conf["detector"]["solid_angle_correction"] = self.solid_angle_correction
        return conf
        
    def set_noise(self, noise=None, noise_spread=None, noise_variation_n=None, noise_filename=None, noise_dataset=None):
        r"""
        Set detector noise type and parameters (this method is called during initialisation)

        Kwargs:
          :noise (str): Noise added to the predicted intensities  (default ``None``)
            
            *Choose one of the following options:*

            ======================= ==================================================================
            ``noise``               Noise model
            ======================= ==================================================================
            ``None``                No noise
            ``'poisson'``           Poisson noise (*shot noise*)
            ``'normal'``            Normal (*Gaussian*) noise 
            ``'uniform'``           Uniformly distributed values within spread limits
            ``'normal_poisson'``    Normal (*Gaussian*) noise on top of Poisson noise (*shot noise*)
            ``'file'``              Noise data from file
            ``'file_poisson'``      Noise data from file on top of Poisson noise (*shot noise*)
            ======================= ==================================================================

          :noise_spread (float): Width (full width at half maximum) of the Gaussian or uniform noise distribution  (default ``None``) 

            .. note:: The argument ``noise_spread`` takes only effect in combination with ``noise='normal'``, ``'uniform'`` or ``'normal_poisson'``

          :noise_filename (str): Location of the HDF5 file that contains the noise data  (default ``None``)

          :noise_dataset (str):  HDF5 dataset (in the file specified by the argument ``noise_filename``) that contains the noise data  (default ``None``)

            .. note:: The arguments ``noise_filename`` and ``noise_dataset`` takes effect only in combination with ``noise='file'`` or ``'file_poisson'``
        """
        if noise in ["file","file_poisson"]:
            self._noise_filename = noise_filename
            self._noise_dataset = noise_dataset
            self._noise = Variation("poisson" if noise == "file_poisson" else None, noise_spread, noise_variation_n, number_of_dimensions=1)
        else:
            self._noise_filename = None
            self._noise_dataset = None
            self._noise = Variation(noise, noise_spread, noise_variation_n, number_of_dimensions=1)

    def set_center_variation(self, center_variation=None, center_spread_x=None, center_spread_y=None, center_variation_n=None):
        """
        Set the variation of the beam center position (this method is called during initialisation)
        
        Kwargs:
          :center_variation(str): Variation of the beam center position (default ``None``)

            *Choose one of the following options:*
            
            ===================== ==============================================
            ``center_variation``  Variation model
            ===================== ==============================================
            ``None``              No variation
            ``'normal'``          Normal (*Gaussian*) random distribution
            ``'uniform'``         Uniform random distribution
            ``'range'``           Equispaced grid around mean center position
            ===================== ==============================================

          :center_spread_x (float): Width of the distribution of center position in *x* [pixel] (default ``None``)
    
          :center_spread_y (float): Width of the distribution of center position in *y* [pixel] (default ``None``)
    
            .. note:: The arguments ``center_spread_y`` and ``center_spread_x`` take effect only in combination with ``center_variation='normal'``, ``'uniform'`` or ``'range'``

          :center_variation_n (int): Number of samples within the specified range (default ``None``)

            .. note:: The argument ``center_variation_n`` takes effect only in combination with ``center_variation='range'``
        """

        
        self._center_variation = Variation(center_variation, [center_spread_x,center_spread_y], center_variation_n, number_of_dimensions=2)
        
    def _init_mask(self, mask, mask_is_cxi_bitmask, mask_filename, mask_dataset, nx, ny, x_gap_size_in_pixel, y_gap_size_in_pixel, cx_hole, cy_hole, hole_diameter_in_pixel):
        if mask is not None or (mask_filename is not None and mask_dataset is not None):
            if mask is not None:
                # Copy mask from array
                self._mask = numpy.array(mask, dtype=numpy.uint16)
            else:
                # Read mask from file
                import h5py
                with h5py.File(mask_filename,"r") as f:
                    self._mask = numpy.array(f[mask_dataset][:,:], dtype=numpy.uint16)
            if not mask_is_cxi_bitmask:
                # Convert maskt to CXI bit format
                self._mask = (self._mask == 0) * PixelMask.PIXEL_IS_MISSING
        elif nx is not None and ny is not None:
            # Initialise empty mask
            self._mask = numpy.zeros(shape=(ny+y_gap_size_in_pixel, nx+x_gap_size_in_pixel),dtype=numpy.uint16)
        else:
            log_and_raise_error(logger, r"Either 'mask' or 'nx' and 'ny' have to be specified.")
            sys.exit(1)
        self._nx = self._mask.shape[1]
        self._ny = self._mask.shape[0]
        # Mask out pixels in gaps
        if y_gap_size_in_pixel > 0:
            cy = int(numpy.ceil((self._ny-1)/2.))
            gy = int(numpy.round(y_gap_size_in_pixel))
            self._mask[cy-gy/2:cy-gy/2+gy,:] |= PixelMask.PIXEL_IS_MISSING
        if x_gap_size_in_pixel > 0:
            cx = int(numpy.ceil((self._nx-1)/2.))
            gx = int(numpy.round(x_gap_size_in_pixel))
            self._mask[:,cx-gx/2:cx-gx/2+gx] |= PixelMask.PIXEL_IS_MISSING
        # Mask out pixels in hole    
        if hole_diameter_in_pixel > 0:
            if cx_hole is None:
                cx_hole = (self._nx-1)/2.
            if cy_hole is None:
                cy_hole = (self._ny-1)/2.
            Y,X = numpy.indices((self._ny,self._nx), dtype=numpy.float64)
            X = X-cx_hole
            Y = Y-cy_hole
            R = numpy.sqrt(X**2 + Y**2)
            tmp = R<=hole_diameter_in_pixel/2.0
            if tmp.sum() > 0:
                self._mask[tmp] |= PixelMask.PIXEL_IS_MISSING

    def get_mask(self,intensities=None, boolmask=False):
        """
        Return mask. The mask has information about the status of each individual detector pixel. The output can be either a CXI bitmask (default) or a boolean mask
    
        For further information and the full bitcode go to :class:`condor.utils.pixelmask.PixelMask`
       
        Kwargs:
          :intensities: Numpy array of photon intensities for masking saturated pixels (default ``None``)

          :boolmask (bool): If ``True`` the output will be a boolean array. Mask values are converted to ``True`` if no bit is set and to ``False`` otherwise
        """
        if intensities is not None:
            if not condor.utils.testing.same_shape(intensities, self._mask):
                log_and_raise_error(logger, "Intensities and mask do not have the same shape")
        M = self._mask.copy()
        if self.saturation_level is not None and intensities is not None:
            M[intensities >= self.saturation_level] |= PixelMask.PIXEL_IS_SATURATED
        if boolmask:
            return numpy.array(M == 0,dtype="bool")
        else:
            return M
    
    def get_cx_mean_value(self):
        """
        Return *x*-coordinate of the mean beam center position
        """
        if self.cx_mean is None:
            return (self._nx-1) / 2.
        else:
            return self.cx_mean

    def get_cy_mean_value(self):
        """
        Return *y*-coordinate of the mean beam center position
        """

        if self.cy_mean is None:
            return (self._ny-1) / 2.
        else:
            return self.cy_mean
        
    def get_next(self):
        """
        Iterate the parameters of the Detector instance and return them as a dictionary
        """
        O = {}
        cx_mean = self.get_cx_mean_value()
        cy_mean = self.get_cy_mean_value()
        cx, cy = self._center_variation.get([cx_mean, cy_mean])
        O["cx"] = cx
        O["cy"] = cy
        O["nx"] = self._nx
        O["ny"] = self._ny
        O["pixel_size"] = self.pixel_size
        O["distance"] = self.distance
        if self.binning is not None:
            O["cx_xxx"] = utils.resample.downsample_pos(cx, self._nx, self.binning)
            O["cy_xxx"] = utils.resample.downsample_pos(cy, self._ny, self.binning)
        return O

    def get_pixel_solid_angle(self, x_off=0., y_off=0.):
        """
        Get the solid angle for a pixel at position ``x_off``, ``y_off`` with respect to the beam center
        
        Kwargs:
          :x_off: *x*-coordinate of the pixel position (center) in unit pixel with respect to the beam center (default 0.)

          :y_off: *y*-coordinate of the pixel position (center) in unit pixel with respect to the beam center (default 0.)
        """
        r_max = numpy.sqrt(x_off**2+y_off**2) * self.pixel_size
        it = isinstance(r_max, collections.Iterable)
        if it:
            r_max = r_max.max()
        if r_max/self.distance < 0.0001:
            # Small angle approximation (fast)
            omega = self.pixel_size**2 / self.distance**2
            if it:
                omega *= numpy.ones_like(r_max)
        else:
            # More precise formula for large angles (slow)
            x_alpha = numpy.arctan2((x_off+0.5)*self.pixel_size, self.distance) - numpy.arctan2((x_off-0.5)*self.pixel_size, self.distance)
            y_alpha = numpy.arctan2((y_off+0.5)*self.pixel_size, self.distance) - numpy.arctan2((y_off-0.5)*self.pixel_size, self.distance)
            omega = 4. * numpy.arcsin(numpy.sin(x_alpha/2.)*numpy.sin(y_alpha/2.))
        return omega

    def get_all_pixel_solid_angles(self, cx, cy):
        """
        Return the solid angles of all detector pixels assuming a beam center at position (``cx``, ``cy``).
        
        Args:
          :cx (float): *x*-coordinate of the center position in unit pixel

          :cy (float): *y*-coordinate of the center position in unit pixel
        """
        Y, X = numpy.meshgrid(numpy.float64(numpy.arange(self._ny))-cy,
                              numpy.float64(numpy.arange(self._nx))-cx,
                              indexing="ij")
        return self.get_pixel_solid_angle(X, Y)
    
    def _get_xy_max_dist(self, cx = None, cy = None, center_variation = False):
        dist_max = []
        for c,dc,n in [[(cx if cx is not None else self.get_cx_mean_value()),self._center_variation.get_spread()[0],self._nx],
                       [(cy if cy is not None else self.get_cy_mean_value()),self._center_variation.get_spread()[1],self._ny]]:
            lim = 0
            if center_variation:
                lim = dc 
                lim = lim*0.5 if lim is not None else 0.
                cv_mode = self._center_variation.get_mode()
                if cv_mode is not None:
                    if "normal" in cv_mode:
                        lim *= 3
            c_min = c - lim/2.
            c_max = c + lim/2.
            res1 = n-1-c_min
            res2 = -c_max
            dist_max.append(res1 if abs(res1) > abs(res2) else res2)
        return dist_max

    def get_p_max_dist(self, cx = None, cy = None, pos = "corner", center_variation = False):
        r"""
        Return 3D position vector of the pixel furthest away from the beam center. If each of the given center position coordinates (``cx``, ``cy``) is ``None`` the beam center is assumed to be located at its mean position
        
        Kwargs:
          :cx (float): *x*-coordinate of the center position in unit pixel (default ``None``)

          :cy (float): *y*-coordinate of the center position in unit pixel (default ``None``)

          :pos (str): Position constraint can be either ``pos='corner'`` or ``pos='edge'``. (default ``'corner'``)

          :center_variation (bool): If ``True`` the beam center variation is taken into account. With respect to the mean position a maximum deviation of *factor/2* times the variational spread is assumed. The *factor* is 3 for Gaussian distributed centers and 1 for others  (default ``False``)
        """
        x, y = self._get_xy_max_dist(cx=cx, cy=cy, center_variation=center_variation)
        xm = x*self.pixel_size
        ym = y*self.pixel_size        
        log_debug(logger, "x = %.1f pix, y = %.1f pix" % (x, y))
        p = numpy.array([0.,0.,self.distance])
        if pos == "corner":
            p[0] = xm
            p[1] = ym
        elif pos == "edge":
            if abs(x) > abs(y):
                p[0] = xm
            else:
                p[1] = ym
        else:
            log_and_raise_error(logger, r"Invalid input: pos=%s. Input must be either 'corner' or 'edge'." % pos)
        return p

    def get_q_max(self, wavelength, cx = None, cy = None, pos = "corner", center_variation = False):
        """
        Return q-vector of maximal lenght

        Args:
          :wavelength (float): Photon wavelength in meters

        Kwargs:
          :cx (float): *x*-coordinate of the center position in unit pixel (default ``None``)

          :cy (float): *y*-coordinate of the center position in unit pixel (default ``None``)
        """
        p = abs(self.get_p_max_dist(cx=cx, cy=cy, pos=pos, center_variation=center_variation))
        q = abs(condor.utils.scattering_vector.q_from_p(p, wavelength))
        return q

    def _get_resolution_element(self, wavelength, cx = None, cy = None, center_variation = False, pos="edge"):
        res = numpy.pi / self.get_q_max(wavelength, cx=cx, cy=cy, pos=pos, center_variation=center_variation)
        return res

    def get_max_resolution(self, wavelength, cx = None, cy = None, center_variation = False):
        """
        Return maximum resolution as a 3D vector (i.e. maximum resolution / momentum transfer in *x*, *y* and *z*)
        
        Args:
          :wavelength (float): Photon wavelength in meters

        Kwargs:
          :cx (float): *x*-coordinate of the center position in unit pixel (default ``None``)

          :cy (float): *y*-coordinate of the center position in unit pixel (default ``None``)

          :center_variation (bool): If ``True`` the beam center variation is taken into account. With respect to the mean position a maximum deviation of *factor/2* times the variational spread is assumed. The *factor* is 3 for Gaussian distributed centers and 1 for others (default ``False``)
        """
        return self._get_resolution_element(wavelength, cx=cx, cy=cy, center_variation=center_variation, pos="corner")

    def get_resolution_element_y(self, wavelength, cx = None, cy = None, center_variation = False):
        """
        Return resolution in *y* in 1/meters

        Args:
          :wavelength (float): Photon wavelength in meters

        Kwargs:
          :cx (float): *x*-coordinate of the center position in unit pixel (default ``None``)

          :cy (float): *y*-coordinate of the center position in unit pixel (default ``None``)

          :center_variation (bool): If ``True`` the beam center variation is taken into account. With respect to the mean position a maximum deviation of *factor/2* times the variational spread is assumed. The *factor* is 3 for Gaussian distributed centers and 1 for others (default ``False``)
        """
        return self._get_resolution_element(wavelength, cx=cx, cy=cy, center_variation=center_variation, pos="corner")[1]

    def get_resolution_element_x(self, wavelength, cx = None, cy = None, center_variation = False):
        """
        Return resolution in *x* in 1/meters

        Args:
          :wavelength (float): Photon wavelength in meters

        Kwargs:
          :cx (float): *x*-coordinate of the center position in unit pixel (default ``None``)

          :cy (float): *y*-coordinate of the center position in unit pixel (default ``None``)

          :center_variation (bool): If ``True`` the beam center variation is taken into account. With respect to the mean position a maximum deviation of *factor/2* times the variational spread is assumed. The *factor* is 3 for Gaussian distributed centers and 1 for others (default ``False``)
        """
        return self._get_resolution_element(wavelength, cx=cx, cy=cy, center_variation=center_variation, pos="corner")[0]

    def get_resolution_element_r(self, wavelength, cx = None, cy = None, center_variation = False):
        """
        Return resolution at the furthes corner position in 1/meters

        Args:
          :wavelength (float): Photon wavelength in meters

        Kwargs:
          :cx (float): *x*-coordinate of the center position in unit pixel (default ``None``)

          :cy (float): *y*-coordinate of the center position in unit pixel (default ``None``)

          :center_variation (bool): If ``True`` the beam center variation is taken into account. With respect to the mean position a maximum deviation of *factor/2* times the variational spread is assumed. The *factor* is 3 for Gaussian distributed centers and 1 for others (default ``False``)
        """
        
        qmax = self.get_q_max(wavelength, cx=cx, cy=cy, pos="corner", center_variation=center_variation)
        res = numpy.pi / length(qmax)
        return res

    def generate_xypix(self, cx=None, cy=None):
        Y, X = numpy.meshgrid(numpy.float64(numpy.arange(self._ny))-(0. if cy is None else cy),
                              numpy.float64(numpy.arange(self._nx))-(0. if cx is None else cx),
                              indexing="ij")
        return X, Y
        
    def generate_qmap(self, wavelength, cx=None, cy=None, extrinsic_rotation=None, order='xyz'):
        X, Y = self.generate_xypix(cx, cy)
        return condor.utils.scattering_vector.generate_qmap(X, Y, self.pixel_size, self.distance, wavelength, extrinsic_rotation=extrinsic_rotation, order=order)

    def generate_qmap_3d(self, wavelength, qn=None, qmax=None, extrinsic_rotation=None, order='xyz'):
        if qn is None and qmax is None:
            qn = max([self._nx, self._ny])
            qmax = self.get_q_max(wavelength, pos="edge")
        elif qn is not None and qmax is not None:
            pass
        else:
            log_and_raise_error(logger, "Either none or both optional arguments qn and qmax have to be passed to this function.")
            return
        return condor.utils.scattering_vector.generate_qmap_3d(qn=qn, qmax=qmax, extrinsic_rotation=extrinsic_rotation, order=order)

    #def generate_rpix_3d(self, qmax, qn, wavelength):
    #    return condor.utils.scattering_vector.generate_rpix_3d(qn, qmax, wavelength, self.distance, self.pixel_size):

    def calculate_polarization_factors(self, cx=None, cy=None, polarization="ignore"):
        if polarization == "ignore":
            P = numpy.ones(shape=(self._ny, self._nx))
        else:
            X, Y = self.generate_xypix(cx=cx, cy=cy)
            P = condor.utils.diffraction.polarization_factor(X, Y, self.distance, polarization=polarization)
        return P
    
    def detect_photons(self, I):
        """
        Return measurement of intensities from an array of expectation values of intensities. This method also returns the mask of the pattern

        Args:
          :I (array): Intensity pattern represented as 2D array
        """
        I_det = self._noise.get(I)
        if self._noise_filename is not None:
            import h5py
            with h5py.File(self._noise_filename,"r") as f:
                ds = f[self._noise_dataset]
                if len(list(ds.shape)) == 2:
                    bg = ds[:,:]
                else:
                    bg = ds[numpy.random.randint(ds.shape[0]),:,:]
            I_det = I_det + bg
        if self.saturation_level is not None:
            I_det = numpy.clip(I_det, -numpy.inf, self.saturation_level)
        if I_det.ndim == 2:
            M_det = self.get_mask(I_det)
        else:
            M_det = None
        return I_det, M_det

    def bin_photons(self, I_det, M_det):
        """
        Return the tuple of binned diffraction pattern and mask. If binning has not been specified a tuple ``(None, None)`` is returned

        Args:
          :I_det (array): Intensity pattern (before binning) represented by a 2D array

          :M_det (array): CXI bitmask (before binning) represented by a 2D array (see also :class:`condor.utils.pixelmask.PixelMask`)
        """
        if self.binning is not None:
            IXxX_det, MXxX_det = utils.resample.downsample(I_det,self.binning,mode="integrate",
                                                           mask2d0=M_det,bad_bits=PixelMask.PIXEL_IS_IN_MASK,min_N_pixels=1)
        else:
            IXxX_det = None
            MXxX_det = None
        return IXxX_det, MXxX_det


        
