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

import sys
import numpy
#from scipy.interpolate import RegularGridInterpolator

import logging
logger = logging.getLogger(__name__)

import condor
import condor.utils.log
from condor.utils.log import log_and_raise_error,log_warning,log_info,log_debug

from condor.utils.variation import Variation
import condor.utils.spheroid_diffraction
import condor.utils.diffraction
import condor.utils.bodies

import condor.utils.emdio

from particle_abstract import AbstractContinuousParticle

ENABLE_MAP_INTERPOLATION = False

class ParticleMap(AbstractContinuousParticle):
    """
    Class for a particle model

    *Model:* Refractive index map sampled on a cubic grid (continuum approximation)

    Args:
      :geometry (str): Geometry type

        *Choose one of the following options:*

          - ``\'custom\'`` - provide map either with an HDF5 file (``map3d_filename``, ``map3d_dataset``) or with a numpy array (``map3d``)

          - ``\'icosahedron\'`` - create map of a uniformly filled icosahedron

          - ``\'cube\'`` - create map of a uniformly filled cube

          - ``\'sphere\'`` - create map of a uniformly filled sphere

          - ``\'spheroid\'`` - create map of a uniformly filled spheroid

      :diameter (float): Particle diameter (not map diameter)

    Kwargs:
      :diameter_variation (str): See :meth:`condor.particle.particle_abstract.AbstractContinuousParticle.set_diameter_variation` (default ``None``)

      :diameter_spread (float): See :meth:`condor.particle.particle_abstract.AbstractContinuousParticle.set_diameter_variation` (default ``None``)

      :diameter_variation_n (int): See :meth:`condor.particle.particle_abstract.AbstractContinuousParticle.set_diameter_variation` (default ``None``)

      :dx: Distance between grid points of the map. This needs to be specified only if ``geometry=`\custom'\'. Depending on whether the geometry is specified by file (``map3d_filename``, ``map3d_dataset``) or by numpy array (``map3d``) for more documentation see :meth:`set_custom_geometry_by_h5file` or :meth:`set_custom_geometry_by_array` respectively (default ``None``)

      :map3d: See :meth:`set_custom_geometry_by_array` (default ``None``)

      :map3d_filename: See :meth:`set_custom_geometry_by_h5file` (default ``None``)

      :map3d_dataset: See :meth:`set_custom_geometry_by_h5file` (default ``None``)

      :emd_id: See :meth:`set_custom_geometry_by_emd_id` (default ``None``)

      :rotation_values (array): See :meth:`condor.particle.particle_abstract.AbstractParticle.set_alignment` (default ``None``)

      :rotation_formalism (str): See :meth:`condor.particle.particle_abstract.AbstractParticle.set_alignment` (default ``None``)

      :rotation_mode (str): See :meth:`condor.particle.particle_abstract.AbstractParticle.set_alignment` (default ``None``)

      :flattening (float): (Mean) value of :math:`a/c`, takes only effect if ``geometry=\'spheroid\'`` (default ``0.75``)
    
      :number (float): Expectation value for the number of particles in the interaction volume. (defaukt ``1.``)

      :arrival (str): Arrival of particles at the interaction volume can be either ``'random'`` or ``'synchronised'``. If ``sync`` at every event the number of particles in the interaction volume equals the rounded value of ``number``. If ``'random'`` the number of particles is Poissonian and ``number`` is the expectation value. (default ``'synchronised'``)

      :position (array): See :class:`condor.particle.particle_abstract.AbstractParticle` (default ``None``)

      :position_variation (str): See :meth:`condor.particle.particle_abstract.AbstractParticle.set_position_variation` (default ``None``)

      :position_spread (float): See :meth:`condor.particle.particle_abstract.AbstractParticle.set_position_variation` (default ``None``)

      :position_variation_n (int): See :meth:`condor.particle.particle_abstract.AbstractParticle.set_position_variation` (default ``None``)

      :material_type (str): See :meth:`condor.particle.particle_abstract.AbstractContinuousParticle.set_material` (default ``\'water\'``)

      :massdensity (float): See :meth:`condor.particle.particle_abstract.AbstractContinuousParticle.set_material` (default ``None``)

      :atomic_composition (dict): See :meth:`condor.particle.particle_abstract.AbstractContinuousParticle.set_material` (default ``None``)          

      :electron_density (float): See :meth:`condor.particle.particle_abstract.AbstractContinuousParticle.set_material` (default ``None``)
    """
    def __init__(self,
                 geometry, diameter = None,
                 diameter_variation = None, diameter_spread = None, diameter_variation_n = None,
                 dx = None,
                 map3d = None, 
                 map3d_filename = None, map3d_dataset = None,
                 emd_id = None,
                 rotation_values = None, rotation_formalism = None, rotation_mode = "extrinsic",
                 flattening = 0.75,
                 number = 1., arrival = "synchronised",
                 position = None, position_variation = None, position_spread = None, position_variation_n = None,
                 material_type = None, massdensity = None, atomic_composition = None, electron_density = None):
        # Initialise base class
        AbstractContinuousParticle.__init__(self,
                                            diameter=diameter, diameter_variation=diameter_variation, diameter_spread=diameter_spread, diameter_variation_n=diameter_variation_n,
                                            rotation_values=rotation_values, rotation_formalism=rotation_formalism, rotation_mode=rotation_mode,                                            
                                            number=number, arrival=arrival,
                                            position=position, position_variation=position_variation, position_spread=position_spread, position_variation_n=position_variation_n,
                                            material_type=material_type, massdensity=massdensity, atomic_composition=atomic_composition, electron_density=electron_density)
        
        # Check for valid geometry
        if geometry not in ["icosahedron", "cube", "sphere", "spheroid", "custom"]:
            log_and_raise_error(logger, "Cannot initialize %s because \'%s\' is not a valid argument for \'geometry\'." % (kwargs["geometry"], self.__class__.__name__))
            sys.exit(1)
        self.geometry = geometry
        
        # Has effect only for spheroids
        self.flattening = flattening

        # Init chache
        self._cache = {}
        self._dx_orig                = None
        self._map3d_orig             = None

        if geometry == "custom":
            if map3d is not None:
                if dx is None:
                    log_and_raise_error(logger, "Cannot initialize custom geometry with \'map3d\' without known grid spacing (\'dx\').")
                    sys.exit(1)
                else:
                    log_debug(logger, "Attempting to initialise custom geometry with \'map3d\'.")
                    if map3d_filename is not None or map3d_dataset is not None or emd_id is not None:
                        log_and_raise_error(logger, "Cannot initialize custom geometry because of ambiguous keyword arguments.")
                        sys.exit(1)
                    self.set_custom_geometry_by_array(map3d, dx)
            elif map3d_filename is not None and map3d_dataset is not None:
                if dx is None:
                    log_and_raise_error(logger, "You are trying to initialise the map with an HDF5 file. You also need to provide the grid spacing \'dx\'")
                    sys.exit(1)
                log_debug(logger, "Attempting to initialise custom geometry with \'map3d_filename\', \'map3d_dataset\' and \'dx\'.")
                if not map3d_filename.endswith(".h5"):
                    log_and_raise_error(logger, "Map file is not an HDF5 file!")
                    sys.exit(1)
                if map3d is not None or emd_id is not None: 
                    log_and_raise_error(logger, "Cannot initialize custom geometry because of ambiguous keyword arguments.")
                    sys.exit(1)
                self.set_custom_geometry_by_h5file(map3d_filename, map3d_dataset, dx)
            elif map3d_filename is not None:
                if not map3d_filename.endswith(".map") and not map3d_filename.endswith(".mrc"):
                    log_and_raise_error(logger, "Map file is not an MRC/MAP file!")
                    sys.exit(1)
                self.set_custom_geometry_by_mrcfile(map3d_filename)
            elif emd_id is not None:
                log_debug(logger, "Attempting to initialise custom geometry with \'emd_id\'.")
                if map3d_filename is not None or map3d_dataset is not None or map3d is not None or dx is not None:
                    log_and_raise_error(logger, "Cannot initialize custom geometry because of ambiguous keyword arguments.")
                    sys.exit(1)
                self.set_custom_geometry_by_emd_id(emd_id)

        if diameter is None:
            self.diameter_mean = self._dx_orig * self._map3d_orig.shape[-1]
            
    def get_conf(self):
        """
        Get configuration in form of a dictionary. Another identically configured ParticleMap instance can be initialised by:

        .. code-block:: python

          conf = P0.get_conf()            # P0: already existing ParticleMap instance
          P1 = condor.ParticleMap(**conf) # P1: new ParticleMap instance with the same configuration as P0  
        """
        conf = {}
        conf.update(AbstractContinuousParticle.get_conf(self))
        conf["geometry"] = self.geometry
        if self.geometry == "custom":
            m,dx = self.get_original_map()
            conf["map3d"] = m
            conf["dx"]    = dx
        if self.geometry == "spheroid":
            conf["flattening"] = self.flattening
        return conf

    def get_next(self):
        """
        Iterate the parameters and return them as a dictionary
        """
        O = AbstractContinuousParticle.get_next(self)
        O["particle_model"] = "map"
        O["geometry"]       = self.geometry
        if self.geometry == "spheroid":
            O["flattening"]     = self.flattening
        return O

    def set_custom_geometry_by_array(self, map3d, dx):
        """
        Set map from numpy array

        Args:
          :map3d (array): 4D numpy array (material index, z, y, x) of float values. If a material is defined (material not ``None``) the values of the map scale the complex refractive index of the material. If no material is defined (materials is ``None``) the map will be casted to complex values and used without any rescaling.

          :dx (float): Grid spacing in unit meter
        """
        # Check shape
        s = numpy.array(map3d.shape)
        if numpy.any(s[-3:]!=s[-1]):
            log_and_raise_error(logger, "Condor only accepts maps with equal spatial dimensions. Current shape is: %s" % str(s[-3:]))
        if self.materials is None:
            # Complex map(s) = refractive index map
            # Check input
            if len(s) == 3:
                map3d = [map3d]
            if len(s) < 3 or len(s) > 4:
                log_and_raise_error(logger, "map3d has %i dimensions but should have 3 or 4." % len(s))
                return
            # Load map(s)
            _map3d = numpy.asarray(map3d)
        else:
            # Real map(s) to be scaled by material's complext refractive index
            # Check input
            if len(s) not in [3,4]:
                log_and_raise_error(logger, "map3d has %i dimensions but it has to have either 3 or 4." % len(s))
                return
            # Load map(s)
            if len(s) == 3:
                n_mat = len(self.materials)
                s = numpy.array([n_mat] + list(s))
                _map3d = numpy.array(n_mat*[map3d], dtype=numpy.float64)
            else:
                if s[0] != len(self.materials):
                    log_and_raise_error(logger, "The first dimension of the map (%i) does not equal the number of specified materials (%i)." % (s[0], len(self.materials)))
                    return
                _map3d = numpy.asarray(map3d, dtype=numpy.float64)
        self._map3d_orig = _map3d
        self._dx_orig    = dx
        self._set_cache(_map3d, dx, geometry="custom")

    def set_custom_geometry_by_h5file(self, map3d_filename, map3d_dataset, dx):
        """
        Load map from dataset in HDF5 file

        If a material is defined (``material_type`` is not ``None``) the absolute values of the map will be rescaled by the complex refractive index of the material. If no material is defined (``material_type=None``) the map will be casted to complex values and used without any rescaling.

        Args:
          :map3d_filename (str): Location of the HDF5 file that contains the map data

          :map3d_dataset (str): Dataset location in the file. The dataset must have three equal dimensions of float values.

          :dx: Grid spacing in unit meter
        """
        import h5py
        with h5py.File(map3d_filename,"r") as f:
            if map3d_dataset is not None:
                ds = map3d_dataset
            elif len(f.keys()) == 1:
                ds = f.keys()[0]
            else:
                log_and_raise_error(logger, "No dataset specified where to find the map.")
            if len(f[ds].shape) == 4:
                map3d = numpy.array(f[ds][:,:,:,:])
            elif len(f[ds].shape) == 3:
                map3d = numpy.array([f[ds][:,:,:]])
            else:
                log_and_raise_error(logger, "Dataset has %i dimensions but it has to have either 3 or 4." % len(f[ds].shape))
                return 
        self.set_custom_geometry_by_array(map3d, dx)                
        
    def set_custom_geometry_by_emd_id(self, emd_id, offset=None, factor=None):
        """
        Fetch map from the EMD by id code.

        The map will be preprocessed by applying an offset and rescaling and by padding the water background with zeros.

        Finally, the avereage value of the map will be rescaled by the refractive index of the associated material.

        Args:
          :emd_id (str): EMD ID code.

          :offset (float): Offset value of the map (MAP = (EM_DATA + OFFSET) X FACTOR)

          :factor (float): Rescale factor of the map (MAP = (EM_DATA + OFFSET) X FACTOR)
        """
        map3d, dx = condor.utils.emdio.fetch_map(emd_id)
        if offset is None and factor is None:
            ed_water = condor.utils.material.AtomDensityMaterial(material_type="water").get_electron_density()
            if len(self.materials) > 1:
                log_and_raise_error(logger, "More than one material defined. This is incompatible with automatic scaling of an EMD map.")
                sys.exit(1)
            ed_particle = self.materials[0].get_electron_density()
            map3d = condor.utils.emdio.preproc_map_auto(map3d, ed_water=ed_water, ed_particle=ed_particle)
        else:
            map3d = condor.utils.emdio.perproc_map_manual(map3d, offset=offset, factor=factor)           
        self.set_custom_geometry_by_array(map3d, dx)                
        
    def set_custom_geometry_by_mrcfile(self, filename, offset=None, factor=None):
        """
        Read map from the MRC file (CCP4 file format, see http://www.ccp4.ac.uk/html/maplib.html).

        The map will be preprocessed by applying an offset and rescaling and by padding the water background with zeros.

        Finally, the avereage value of the map will be rescaled by the refractive index of the associated material.

        Args:
          :filename (str): Filename of MRC file.

          :offset (float): Offset value of the map (MAP = (EM_DATA + OFFSET) X FACTOR)

          :factor (float): Rescale factor of the map (MAP = (EM_DATA + OFFSET) X FACTOR)
        """
        map3d, dx = condor.utils.emdio.read_map(filename)
        if offset is None and factor is None:
            ed_water = condor.utils.material.AtomDensityMaterial(material_type="water").get_electron_density()
            if len(self.materials) > 1:
                log_and_raise_error(logger, "More than one material defined. This is incompatible with automatic scaling of an EMD map.")
                sys.exit(1)
            ed_particle = self.materials[0].get_electron_density()
            map3d = condor.utils.emdio.preproc_map_auto(map3d, ed_water=ed_water, ed_particle=ed_particle)
        else:
            map3d = condor.utils.emdio.perproc_map_manual(map3d, offset=offset, factor=factor) 
        self.set_custom_geometry_by_array(map3d, dx)                

        
    def get_new_dn_map(self, O, dx_required, dx_suggested, photon_wavelength):
        """
        Return the a new refractive index map

        Args:

          :O (dict): Parameter dictionary as returned from :meth:`condor.particle.particle_map.get_next`

          :dx_required (float): Required resolution (grid spacing) of the map. An error is raised if the resolution of the map has too low resolution

          :dx_suggested (float): Suggested resolution (grid spacing) of the map. If the map has a very high resolution it will be interpolated to a the suggested resolution value

          :photon_wavelength (float): Photon wavelength in unit meter 
        """
        m,dx = self.get_new_map(O=O, dx_required=dx_required, dx_suggested=dx_suggested)
        if self.materials is not None:
            dn = numpy.zeros(shape=(m.shape[1], m.shape[2], m.shape[3]), dtype=numpy.complex128)
            for mat_i, m_i in zip(self.materials, m):
                dn_i = mat_i.get_dn(photon_wavelength=photon_wavelength)
                dn += m_i * dn_i
        else:
            dn = m[0]    
        return dn,dx

    def get_current_map(self):
        """
        Return the current map
        """
        return self._map3d, self._dx

    def get_original_map(self):
        """
        Return the original map

        """
        return self._map3d_orig, self._dx_orig
    
    def _get_map3d(self, O = None, dx_required = None, dx_suggested = None): 
        if O is not None:
            m,dx = self.get_new_map(O, dx_required, dx_suggested)
        else:
            m,dx = self.get_current_map()
        return m, dx

    def _set_cache(self, map3d, dx, geometry, diameter=None, flattening=None):
        self._cache = {
            "map3d"      : map3d,
            "dx"         : dx,
            "geometry"   : geometry,
            "diameter"   : diameter,
            "flattening" : flattening,
        }
    
    def _is_map_in_cache(self, O, dx_required):
        # Empty cache?
        if not self._cache:
            return False
        # Custom map?
        elif O["geometry"] == "custom":
            return False
        # Correct geometry?
        elif self._cache["geometry"] != O["geometry"]:
            return False
        # Correct size?
        elif abs(self._cache["diameter"] - O["diameter"]) > 1E-10:
            return False
        # Correct spheroid flattening?
        elif O["geometry"] == "spheroid":
            if abs(self._cache["flattening"] - O["flattening"]) > 1E-10:
                return False
        # Sufficient resolution?
        elif self._cache["dx"] > dx_required:
            return False
        else:
            return True
        
    def get_new_map(self, O, dx_required, dx_suggested):
        """
        Return new map with given parameters

        Args:

          :O (dict): Parameter dictionary as returned from :meth:`condor.particle.particle_map.get_next`

          :dx_required (float): Required resolution (grid spacing) of the map. An error is raised if the resolution of the map has too low resolution

          :dx_suggested (float): Suggested resolution (grid spacing) of the map. If the map has a very high resolution it will be interpolated to a the suggested resolution value
        """
        
        if O["geometry"] in ["icosahedron", "sphere", "spheroid", "cube"]:
            
            if not self._is_map_in_cache(O, dx_required):

                dx = dx_suggested
                n_mat = len(self.materials)
                
                if O["geometry"] == "icosahedron":
                    m_tmp = self._get_map_icosahedron(O["diameter"]/2., dx)

                elif O["geometry"] == "spheroid":
                    a = condor.utils.spheroid_diffraction.to_spheroid_semi_diameter_a(O["diameter"],O["flattening"])
                    c = condor.utils.spheroid_diffraction.to_spheroid_semi_diameter_c(O["diameter"],O["flattening"])
                    m_tmp = self._get_map_spheroid(a, c, dx)

                elif O["geometry"] == "sphere":
                    m_tmp = self._get_map_sphere(O["diameter"]/2., dx)

                elif O["geometry"] == "cube":
                    m_tmp = self._get_map_cube(O["diameter"], dx)

                else:
                    log_and_raise_error(logger, "Particle map geometry \"%s\" is not implemented. Change your configuration and try again." % O["geometry"])
                    sys.exit(1)

                m = numpy.array(n_mat * [m_tmp])
                    
                self._set_cache(map3d=m,
                                dx=dx,
                                geometry=O["geometry"],
                                diameter=O["diameter"],
                                flattening=(None if O["geometry"] != "spheroid" else O["flattening"]))

            else:

                log_debug(logger, "No need for calculating a new map. Reading map from cache.")
                m  = self._cache["map3d"]
                dx = self._cache["dx"]

        elif O["geometry"] == "custom":

            rescale_factor = O["diameter"] / self.diameter_mean
            dx_rescaled = self._cache["dx"] * rescale_factor
            
            # Current map too coarsely sampled?
            if (dx_rescaled/dx_required > 1.) and not numpy.isclose(dx_rescaled/dx_required, 1.):

                # Cached map (original) also too coarsely sampled? 
                if self._dx_orig/dx_required > 1. and not numpy.isclose(self._dx_orig/dx_required, 1.):
                    # Not fine enough -> exit
                    log_and_raise_error(logger, "Resolution of given custom map is insufficient for simulation. Required is at most %e m vs. provided %e m." % (dx_required, self._dx_orig))
                    sys.exit(1)
                else:
                    # Change back to original fine map
                    self._set_cache(map3d=self._map3d_orig,
                                    dx=self._dx_orig,
                                    geometry="custom")
                    
            # Can we downsample current map?
            # MAX: We would do this only for performance reasons but have not found a good way of downsampling without introducing artifacts
            #if (dx_suggested/dx_rescaled >= 2.) and (dx_suggested/self._dx_orig >= 2.) and ENABLE_MAP_INTERPOLATION:
            #    print "ENABLE_MAP_INTERPOLATION=%i" % ENABLE_MAP_INTERPOLATION
            #    N1 = self._map3d_orig.shape[0]
            #    m1 = numpy.zeros(shape=(N1,N1,N1), dtype=numpy.float64)
            #    m1[:self._map3d_orig.shape[0],:self._map3d_orig.shape[0],:self._map3d_orig.shape[0]] = self._map3d_orig[:,:,:]
            #    fm1 = numpy.fft.fftshift(numpy.fft.ifftn(m1))
            #    N1 = m1.shape[0]
            #    N2 = int(numpy.ceil(N1 * self._dx_orig / dx_suggested))
            #    x1 = numpy.linspace(-0.5,0.5,N2)*(1-0.5/N2)
            #    Z,Y,X = numpy.meshgrid(x1,x1,x1,indexing="ij")
            #    coords = numpy.array([[z,y,x] for z,y,x in zip(Z.ravel(),Y.ravel(),X.ravel())])
            #    m2 = abs(numpy.fft.fftshift(condor.utils.nfft.nfft(fm1,coords).reshape((N2,N2,N2))))
            #    #from pylab import *
            #    #imsave("m1.png", m1.sum(0))
            #    #imsave("m2.png", m2.sum(0))
            #    self._dx    = self._dx_orig * float(N1)/float(N2)
            #    self._map3d = m2 / m2.sum() * m1.sum()

            m  = self._cache["map3d"]
            dx = rescale_factor * self._cache["dx"]
            
        return m,dx
                
    def _get_map_sphere(self, radius, dx):
        nR = radius/dx
        N = int(round((nR*1.2)*2))
        m = condor.utils.bodies.make_sphere_map(N,nR)
        return numpy.asarray(m, dtype=numpy.float64)
 
    def _get_map_spheroid(self, a, c, dx, rotation=None):
        # maximum radius
        Rmax = max([a,c])
        # maximum radius in pixel
        nRmax = Rmax/dx
        # dimensions in pixel
        nA = a/dx
        nC = c/dx
        # leaving a bit of free space around spheroid
        N = int(round((nRmax*1.2)*2))
        m = condor.utils.bodies.make_spheroid_map(N,nA,nC,rotation)
        return numpy.asarray(m, dtype=numpy.float64)

    def _get_map_icosahedron(self, radius, dx):
        # icosahedon size parameter
        a = radius*(16*numpy.pi/5.0/(3+numpy.sqrt(5)))**(1/3.0)
        # radius at corners in meter
        Rmax = numpy.sqrt(10.0+2*numpy.sqrt(5))*a/4.0 
        # radius at corners in pixel
        nRmax = Rmax/dx 
        # leaving a bit of free space around icosahedron 
        N = int(numpy.ceil(2.3*(nRmax)))
        log_info(logger,"Building icosahedron with radius %e (%i pixel) in %i x %i x %i voxel cube." % (radius,nRmax,N,N,N))
        m = condor.utils.bodies.make_icosahedron_map(N,nRmax)
        return numpy.asarray(m, dtype=numpy.float64)
        
    def _get_map_cube(self, a, dx):
        # edge_length in pixels
        nel = a/dx 
        # leaving a bit of free space around
        N = int(numpy.ceil(2.3*nel))
        # make map
        X,Y,Z = 1.0*numpy.mgrid[0:N,0:N,0:N]
        X = X - (N-1)/2.
        Y = Y - (N-1)/2.
        Z = Z - (N-1)/2.
        DX = abs(X)-nel/2.
        DY = abs(Y)-nel/2.
        DZ = abs(Z)-nel/2.
        D = numpy.array([DZ,DY,DX])
        m = numpy.zeros(shape=(N,N,N))
        Dmax = D.max(0)
        m[Dmax<-0.5] = 1.
        temp = (Dmax<0.5)*(m==0.)
        d = Dmax[temp]
        m[temp] = 0.5-d
        return numpy.asarray(m, dtype=numpy.float64)
