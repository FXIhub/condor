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

import sys
import numpy
import scipy.stats
from scipy.interpolate import RegularGridInterpolator

import logging
logger = logging.getLogger(__name__)

import condor
import condor.utils.log
from condor.utils.log import log_and_raise_error,log_warning,log_info,log_debug

from condor.utils.variation import Variation
import condor.utils.spheroid_diffraction
import condor.utils.diffraction
import condor.utils.bodies

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

      :rotation_values (array): See :meth:`condor.particle.particle_abstract.AbstractParticle.set_alignment` (default ``None``)

      :rotation_formalism (str): See :meth:`condor.particle.particle_abstract.AbstractParticle.set_alignment` (default ``None``)

      :rotation_mode (str): See :meth:`condor.particle.particle_abstract.AbstractParticle.set_alignment` (default ``None``)

      :flattening (float): (Mean) value of :math:`a/c`, takes only effect if ``geometry=\'spheroid\'`` (default ``0.75``)
    
:number_density (float): Number density of this particle species in units of the interaction volume. (defaukt ``1.``)

      :arrival (str): Arrival of particles at the interaction volume can be either ``'random'`` or ``'synchronised'``. If ``sync`` at every event the number of particles in the interaction volume equals the rounded value of the ``number_density``. If ``'random'`` the number of particles is Poissonian and the ``number_density`` is the expectation value. (default ``'synchronised'``)

      :position (array): See :class:`condor.particle.particle_abstract.AbstractParticle` (default ``None``)

      :position_variation (str): See :meth:`condor.particle.particle_abstract.AbstractParticle.set_position_variation` (default ``None``)

      :position_spread (float): See :meth:`condor.particle.particle_abstract.AbstractParticle.set_position_variation` (default ``None``)

      :position_variation_n (int): See :meth:`condor.particle.particle_abstract.AbstractParticle.set_position_variation` (default ``None``)

      :material_type (str): See :meth:`condor.particle.particle_abstract.AbstractContinuousParticle.set_material` (default ``\'water\'``)

      :massdensity (float): See :meth:`condor.particle.particle_abstract.AbstractContinuousParticle.set_material` (default ``None``)

      :atomic_composition (dict): See :meth:`condor.particle.particle_abstract.AbstractContinuousParticle.set_material` (default ``None``)          
    """
    def __init__(self,
                 geometry, diameter,
                 diameter_variation = None, diameter_spread = None, diameter_variation_n = None,
                 dx = None,
                 map3d = None, 
                 map3d_filename = None, map3d_dataset = None,
                 rotation_values = None, rotation_formalism = None, rotation_mode = "extrinsic",
                 flattening = 0.75,
                 number_density = 1., arrival = "synchronised",
                 position = None, position_variation = None, position_spread = None, position_variation_n = None,
                 material_type = None, massdensity = None, atomic_composition = None):
        # Initialise base class
        AbstractContinuousParticle.__init__(self,
                                            diameter=diameter, diameter_variation=diameter_variation, diameter_spread=diameter_spread, diameter_variation_n=diameter_variation_n,
                                            rotation_values=rotation_values, rotation_formalism=rotation_formalism, rotation_mode=rotation_mode,                                            
                                            number_density=number_density, arrival=arrival,
                                            position=position, position_variation=position_variation, position_spread=position_spread, position_variation_n=position_variation_n,
                                            material_type=material_type, massdensity=massdensity, atomic_composition=atomic_composition)
        
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
            if (map3d_filename is not None or map3d_dataset is not None) and map3d is not None:
                log_and_raise_error(logger, "Cannot initialize custom geometry because the keyword arguments \'map3d_filename\' or \'map3d_dataset\' and \'map3d\' are not None.")
                sys.exit(1)
            elif (map3d_filename is None or map3d_dataset is None) and map3d is None:
                log_and_raise_error(logger, "Cannot initialize custom geometry because the keyword arguments for either (\'map3d_filename\' or \'map3d_dataset\') or \'map3d\' are None.")
                sys.exit(1)
            elif dx is None:
                log_and_raise_error(logger, "Cannot initialize custom geometry because the keyword argument \'dx\' is None.")
                sys.exit(1)                
            elif map3d is not None:
                self.set_custom_geometry_by_array(map3d, dx)
            elif (map3d_filename is not None and map3d_dataset is not None):
                self.set_custom_geometry_by_h5file(map3d_filename, map3d_dataset, dx)
            
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
          :map3d (array): Numpy array with three equal dimensions of float values. If a material is defined (``material_type`` is not ``None``) the absolute values of the map will be rescaled by the complex refractive index of the material. If no material is defined (``material_type=None``) the map will be casted to complex values and used without any rescaling.

          :dx (float): Grid spacing in unit meter
        """
        s = numpy.array(map3d.shape)
        if not numpy.all(s==s[0]):
            log_and_raise_error(logger, "Condor only accepts maps with equal dimensions.")
            return
        self._map3d_orig = map3d
        self._dx_orig    = dx
        self._map3d      = map3d
        self._set_cache(map3d, dx, geometry="custom")

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
            map3d = numpy.array(f[ds][:,:,:])
        self.set_custom_geometry_by_array(map3d, dx)                

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
    
    def _is_map_in_cache(self, O):
        # Empty cache?
        if not self._cache:
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
            
            if not self._is_map_in_cache(O):

                dx = dx_suggested
                
                if O["geometry"] == "icosahedron":
                    m = self._get_map_icosahedron(O["diameter"]/2., dx)

                elif O["geometry"] == "spheroid":
                    a = condor.utils.spheroid_diffraction.to_spheroid_semi_diameter_a(O["diameter"],O["flattening"])
                    c = condor.utils.spheroid_diffraction.to_spheroid_semi_diameter_c(O["diameter"],O["flattening"])
                    m = self._get_map_spheroid(a, c, dx)

                elif O["geometry"] == "sphere":
                    m = self._put_sphere(O["diameter"]/2., dx)

                elif O["geometry"] == "cube":
                    m = self._get_map_cube(O["diameter"]/2., dx)

                else:
                    log_and_raise_error(logger, "Particle map geometry \"%s\" is not implemented. Change your configuration and try again." % O["geometry"])
                    sys.exit(1)

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
            
            dx_needed = O["diameter"] / self.diameter_mean * self._cache["dx"]
            
            # Map fine enough?
            if dx_needed/dx_required >= 1.:
            
                if self._dx_orig/dx_required >= 1.:
                    # Not fine enough -> exit
                    log_and_raise_error(logger, "Resolution of given custom map is insufficient for simulation. required %e m vs. provided %e m." % (dx_required, self._dx_orig))
                    sys.exit(1)
                    
                # Change back to original fine map
                self.set_cache(map3d=self._map3d_orig,
                               dx=self._dx_orig,
                               geometry="custom")
                    
            # Can we downsample current map?
            #if (dx_suggested/dx_needed >= 2.) and (dx_suggested/self._dx_orig >= 2.) and ENABLE_MAP_INTERPOLATION:
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
            dx = O["diameter"] / self.diameter_mean * self._cache["dx"]
            
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
