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

ENABLE_MAP_INTERPOLATION = True

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
        self._geometry               = None
        self._dx                     = None
        self._dx_orig                = None
        self._map3d                  = None
        self._map3d_orig             = None
        self._diameter               = None
        self._flattening             = None

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
        self._dx_orig = dx
        self._map3d = map3d
        self._dx = dx
        
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
            map3d = numpy.array(f[ds][:,:,:], dtype="float")
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
            
    def get_new_map(self, O, dx_required, dx_suggested):
        """
        Return new map with given parameters

        Args:

          :O (dict): Parameter dictionary as returned from :meth:`condor.particle.particle_map.get_next`

          :dx_required (float): Required resolution (grid spacing) of the map. An error is raised if the resolution of the map has too low resolution

          :dx_suggested (float): Suggested resolution (grid spacing) of the map. If the map has a very high resolution it will be interpolated to a the suggested resolution value
        """
        if O["geometry"] != "custom":
            # Decide whether we need to build a new map
            build_map = False
            if self._map3d is None:
                build_map = True
            if self._diameter is None:
                build_map = True
            else:
                if abs(self._diameter - O["diameter"]) > 1E-10:
                    build_map = True
                if O["flattening"] is not None and self._flattening is not None:
                    if abs(self._flattening - O["flattening"]) > 1E-10:
                        build_map = True
            if self._dx > dx_required:
                build_map = True
                self._map3d = None
                
            if build_map:
                self._map3d = None
                self._dx = dx_suggested
                if O["geometry"] == "icosahedron":
                    self._put_icosahedron(O["diameter"]/2.)
                elif O["geometry"] == "spheroid":
                    a = condor.utils.spheroid_diffraction.to_spheroid_semi_diameter_a(O["diameter"],O["flattening"])
                    c = condor.utils.spheroid_diffraction.to_spheroid_semi_diameter_c(O["diameter"],O["flattening"])
                    self._put_spheroid(a, c)
                elif O["geometry"] == "sphere":
                    self._put_sphere(O["diameter"]/2.)
                elif O["geometry"] == "cube":
                    self._put_cube(O["diameter"]/2.)
                else:
                    log_and_raise_error(logger, "Particle map geometry \"%s\" is not implemented. Change your configuration and try again." % O["geometry"])
                    sys.exit(1)
            m = self._map3d
            dx = self._dx
            
        else:
            dx = O["diameter"] / self.diameter_mean * self._dx
            # Map fine enough?
            if dx/dx_required >= 1.:
                # Original map fine enough?
                if self._dx_orig/dx_required >= 1.:
                    # Not fine enough -> exit
                    log_and_raise_error(logger, "Resolution of custom map is not sufficient for simulation. %e, %e" % (dx_required, self._dx))
                    sys.exit(1)
                else:
                    # Change back to original fine map
                    self._map3d = self._map3d_orig
            # Can we downsample current map?
            if (dx_suggested/dx >= 2.) and (dx_suggested/self._dx_orig >= 2.) and ENABLE_MAP_INTERPOLATION:
                N1 = self._map3d_orig.shape[0]
                m1 = numpy.zeros(shape=(N1,N1,N1), dtype=numpy.float64)
                m1[:self._map3d_orig.shape[0],:self._map3d_orig.shape[0],:self._map3d_orig.shape[0]] = self._map3d_orig[:,:,:]
                fm1 = numpy.fft.fftshift(numpy.fft.ifftn(m1))
                N1 = m1.shape[0]
                N2 = int(numpy.ceil(N1 * self._dx_orig / dx_suggested))
                x1 = numpy.linspace(-0.5,0.5,N2)*(1-0.5/N2)
                Z,Y,X = numpy.meshgrid(x1,x1,x1,indexing="ij")
                coords = numpy.array([[z,y,x] for z,y,x in zip(Z.ravel(),Y.ravel(),X.ravel())])
                m2 = abs(numpy.fft.fftshift(condor.utils.nfft.nfft(fm1,coords).reshape((N2,N2,N2))))
                #from pylab import *
                #imsave("m1.png", m1.sum(0))
                #imsave("m2.png", m2.sum(0))
                self._dx    = self._dx_orig * float(N1)/float(N2)
                self._map3d = m2/m2.sum()*m1.sum()
            dx = O["diameter"] / self.diameter_mean * self._dx
            m = self._map3d

        self._geometry   = O["geometry"]
        self._diameter   = O["diameter"]
        self._flattening = O.get("flattening",None)
        return m,dx
            
    def _put_custom_map(self, map_add):
        self._map3d = numpy.array(map_add, dtype=numpy.float64)
    
    def _put_sphere(self, radius):
        nR = radius/self._dx
        N = int(round((nR*1.2)*2))
        spheremap = condor.utils.bodies.make_sphere_map(N,nR)
        self._put_custom_map(spheremap)
 
    def _put_spheroid(self, a, c, rotation=None):
        # maximum radius
        Rmax = max([a,c])
        # maximum radius in pixel
        nRmax = Rmax/self._dx
        # dimensions in pixel
        nA = a/self._dx
        nC = c/self._dx
        # leaving a bit of free space around spheroid
        N = int(round((nRmax*1.2)*2))
        spheromap = condor.utils.bodies.make_spheroid_map(N,nA,nC,rotation)
        self._put_custom_map(spheromap)

    def _put_icosahedron(self, radius):
        # icosahedon size parameter
        a = radius*(16*numpy.pi/5.0/(3+numpy.sqrt(5)))**(1/3.0)
        # radius at corners in meter
        Rmax = numpy.sqrt(10.0+2*numpy.sqrt(5))*a/4.0 
        # radius at corners in pixel
        nRmax = Rmax/self._dx 
        # leaving a bit of free space around icosahedron 
        N = int(numpy.ceil(2.3*(nRmax)))
        log_info(logger,"Building icosahedron with radius %e (%i pixel) in %i x %i x %i voxel cube." % (radius,nRmax,N,N,N))
        icomap = condor.utils.bodies.make_icosahedron_map(N,nRmax)
        self._put_custom_map(icomap)

    def _put_cube(self, a):
        # edge_length in pixels
        nel = a/self._dx 
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
        cubemap = numpy.zeros(shape=(N,N,N))
        Dmax = D.max(0)
        cubemap[Dmax<-0.5] = 1.
        temp = (Dmax<0.5)*(cubemap==0.)
        d = Dmax[temp]
        cubemap[temp] = 0.5-d
        self._put_custom_map(cubemap)        

