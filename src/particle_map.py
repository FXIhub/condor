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

import logging
logger = logging.getLogger("Condor")
import utils.log
from utils.log import log 

from utils.variation import Variation
import utils.spheroid_diffraction
import utils.diffraction
import utils.bodies

from particle_abstract import AbstractContinuousParticleSpecies

class ParticleSpeciesMap(AbstractContinuousParticleSpecies):
    def __init__(self,**kwargs):
        """
        Function initializes ParticleSpeciesMap object:
        =============================================

        Arguments:

        Keyword arguments (if not given variable \'X\' is set to default value \'[X_default]\'):
        
        EITHER (default):
        - N_fine: Edge length in number of voxels [1]

        OR:
        - map3d_fine: Cubic 3d map.
          EITHER:

          - dx_fine: Real space sampling distance.
          OR:
          - oversampling_fine: Real space oversampling in relation to detector/source configuration. [1]
        
        OR:
        - geometry: Geometry that shall be generated (icosahedron, sphere, spheroid, cube)
          - oversampling_fine: Additional oversampling for the initial map self.map3d_fine [1.]
          ADDITIONAL (spheroid):
          - flattening

        ADDITIONAL:
        - euler_angle_0, euler_angle_1, euler_angle_2: Euler angles defining orientation of 3D grid in the experimental reference frame (beam axis). [0.0,0.0,0.0]
        - size: Characteristic size of the object in meters. [1]
        - parent: Input object that this SampleMap object shall be linked to. [None]
        - SAMPLING:
          EITHER:
          - dX_fine: Real space sampling distance.
          OR:
          - oversampling_fine: Real space oversampling in relation to detector/source configuration. [1]
        - MATERIAL:
          (For more documentation look at the help of the Material class.)
          - material_type: predefined material type. 
          OR
          - massdensity: massdensity of the component
          - cX, cY, ... : atomic composition

        """
        # Check for valid set of keyword arguments
        self.req_keys += ["geometry"]
        self.opt_keys += ["flattening","flattening_variation","flattening_spread","flattening_variation_n"]
        if kwargs["geometry"] is not in ["icosahedron", "cube", "sphere", "spheroid", "custom"]:
            log(logger.error,"Cannot initialize %s because \'%s\' is not a valid argument for \'geometry\'." % (kwargs["geometry"], self.__class__.__name__))
            sys.exit(1)
        # Start initialisation
        AbstractContinuousParticleSpecies.__init__(self,**kwargs)
        self.geometry = kwargs["geometry"]
        self.flattening_mean = kwargs.get("flattening",1.)
        self.set_flattening_variation(flattening_variation=kwargs.get("flattening_variation",None),flattening_spread=kwargs.get("flattening_spread",None),flattening_variation_n=kwargs.get("flattening_variation_n",None))
        self.geometry = kwargs["geometry"]

        # Init chache
        self._old_map3d_diameter               = None
        self._old_map3d_geometry               = None
        self._old_map3d_dx                     = None
        self._old_map3d                        = None

    def get_next(self):
        O = AbstractContinuousParticleSpecies.get_next(self)
        O["flattening"] = self._get_next_flattening()
        return O

    def set_flattening_variation(self,flattening_variation=None,flattening_spread=None,flattening_variation_n=None,**kwargs):
        self._flattening_variation = Variation(flattening_variation,flattening_spread,flattening_variation_n,name="spheroid flattening")       

    def _get_next_flattening(self):
        f = self._flattening_variation.get(self.flattening_mean)
        # Non-random 
        if self._flattening_variation._mode in [None,"range"]:
            if f <= 0:
                log(logger.error,"Spheroid flattening smaller-equals zero. Change your configuration.")
            else:
                return f
        # Random 
        else:
            if f <= 0.:
                log(logger.warning,"Spheroid flattening smaller-equals zero. Try again.")
                return self._get_next_flattening()
            else:
                return f

    def get_map3d(self,O,dx_required,dx_suggested):
        self._build_map(O,dx_required,dx_suggested)
        dx = self._old_map3d_dx
        return m, dx
            
    def _build_map(self,O,dx_required,dx_suggested):
        build_map = False
        if self._old_map3d is None:
            build_map = True
        if self._old_map3d_diameter is None:
            build_map = True
        else:
            if abs(self._old_map3d_diameter - O["diameter"]) > 1E-10:
                build_map = True
        if self._old_map3d_dx > dx_required:
            build_map = True
            self._old_map3d = None
        if not build_map:
            return

        self._old_map3d = None
        self._old_map3d_dx = dx_suggested
        self._old_map3d_diameter = O["diameter"]
        
        if O["geometry"] is not "custom":
            n = self.material.get_dn()
            if O["geometry"] == "icosahedron":
                self._put_icosahedron(O["diameter"]/2., dn)
            elif O["geometry"] == "spheroid":
                a = utils.spheroid_diffraction.to_spheroid_semi_diameter_a(O["diameter"],O["flattening"])
                c = utils.spheroid_diffraction.to_spheroid_semi_diameter_c(O["diameter"],O["flattening"])
                self._put_spheroid(a, c, n)
            elif O["geometry"] == "sphere":
                self._put_sphere(O["diameter"]/2., dn)
            elif O["geometry"] == "cube":
                self._put_cube(O["diameter"]/2., dn)
            else:
                log(logger.error,"Particle map geometry \"%s\" is not implemented. Change your configuration and try again." % O["geometry"])
                sys.exit(1)
        else:
            if "map3d" in O:
                s = numpy.array(O["map3d"].shape)
                if not numpy.all(s==s[0]):
                    log(logger.error,"Condor only accepts maps with equal dimensions.")
                    return
                self._old_map3d = O['map3d']
                self._old_map3d_dx = O["diameter"]/float(s[0]-1)
            elif "filename" in O:
                import h5py
                with h5py.File(O["filename"],"r") as f:
                    if len(f.items()) == 1:
                        d = f.items()[0][1][:,:,:]
                    else:
                        d = f["data"][:,:,:]
                    s = numpy.array(d.shape)
                    if not numpy.all(s==s[0]):
                        log(logger.error,"Condor only accepts maps with equal dimensions.")
                        sys.exit(0)
                    self._old_map3d = d
                    self._old_map3d_dx = O["diameter"]/float(s[0]-1)
            else:
                log(logger.error,"For a custom geometry either map3d or filename has to be specified to load the map. Change your configuration and try again.")
                sys.exit(1)

    #def _put_custom_map(self,map_add,**kwargs):
    #    unit = kwargs.get("unit","meter")
    #    p = numpy.array([kwargs.get("z",0.),kwargs.get("y",0.),kwargs.get("x",0.)])
    #    origin = kwargs.get("origin","middle")
    #    mode = kwargs.get("mode","factor")
    #    dn = kwargs.get("n",None)
    #    if dn == None:
    #        factor = 1.
    #    else:
    #        factor = abs(dn)/abs(self._get_dn())
    #    if self._old_map3d == None:
    #        self._old_map3d = numpy.array(map_add,dtype="float64")
    #    else:
    #        self._old_map3d = utils.bodies.array_to_array(map_add,self._old_map3d,p,origin,mode,0.,factor)

    def _put_custom_map(self, map_add, n):
        self._old_map3d = numpy.array(map_add, dtype=numpy.complex128) * n
    
    def _put_sphere(self, radius, n):
        nR = radius/self._old_map3d_dx
        N = int(round((nR*1.2)*2))
        spheremap = utils.bodies.make_sphere_map(N,nR)
        self._put_custom_map(spheremap, n)
 
    def _put_spheroid(self, a, c, n, e0=0., e1=0., e2=0.):
        # maximum radius
        Rmax = max([a,c])
        # maximum radius in pixel
        nRmax = Rmax/self._old_map3d_dx
        # dimensions in pixel
        nA = a/self._old_map3d_dx
        nC = c/self._old_map3d_dx
        # leaving a bit of free space around spheroid
        N = int(round((nRmax*1.2)*2))
        spheromap = utils.bodies.make_spheroid_map(N,nA,nC,e0,e1,e2)
        self._put_custom_map(spheromap, n)

    def _put_icosahedron(self, radius, n, e0=0., e1=0., e2=0.):
        # icosahedon size parameter
        a = radius*(16*numpy.pi/5.0/(3+numpy.sqrt(5)))**(1/3.0)
        # radius at corners in meter
        Rmax = numpy.sqrt(10.0+2*numpy.sqrt(5))*a/4.0 
        # radius at corners in pixel
        nRmax = Rmax/self._old_map3d_dx 
        # leaving a bit of free space around icosahedron 
        N = int(numpy.ceil(2.3*(nRmax)))
        log(logger.info,"Building icosahedron with radius %e (%i pixel) in %i x %i x %i voxel cube." % (radius,nRmax,N,N,N))
        icomap = utils.bodies.make_icosahedron_map(N,nRmax,e0,e1,e2)
        self._put_custom_map(icomap, n)

    def _put_cube(self, a, n, e0=0., e1=0., e2=0.):
        # edge_length in pixels
        nel = a/self._old_map3d_dx 
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
        self._put_custom_map(cubemap, n)        

    def plot_map3d(self,mode='surface'):
        try:
            from enthought.mayavi import mlab
        except:
            from mayavi import mlab
        if mode=='planes':
            s = mlab.pipeline.scalar_field(abs(self._old_map3d))
            plane_1 = mlab.pipeline.image_plane_widget(s,plane_orientation='x_axes',
                                                       slice_index=self._old_map3d.shape[2]/2)
            plane_2 = mlab.pipeline.image_plane_widget(s,plane_orientation='y_axes',
                                                       slice_index=self._old_map3d.shape[1]/2)
            mlab.show()
        elif mode=='surface':
            mlab.contour3d(abs(self._old_map3d))
        else:
            log(logger.error,"No valid mode given.")
            
    def plot_fmap3d(self):
        from enthought.mayavi import mlab
        import spimage
        M = spimage.sp_image_alloc(self._old_map3d.shape[0],self._old_map3d.shape[1],self._old_map3d.shape[2])
        M.image[:,:,:] = self._old_map3d[:,:,:]
        M.mask[:,:,:] = 1
        fM = spimage.sp_image_fftw3(M)
        fM.mask[:,:,:] = 1
        fsM = spimage.sp_image_shift(fM)
        self.fmap3d = abs(fsM.image).copy()
        self.fmap3d[self.fmap3d!=0] = numpy.log10(self.fmap3d[self.fmap3d!=0])
        
        s = mlab.pipeline.scalar_field(self.fmap3d)
        plane_1 = mlab.pipeline.image_plane_widget(s,plane_orientation='x_axes',
                                                   slice_index=self.fmap3d.shape[2]/2)
        plane_2 = mlab.pipeline.image_plane_widget(s,plane_orientation='y_axes',
                                                   slice_index=self.fmap3d.shape[1]/2)
        mlab.show()


    def save_map3d(self,filename):
        """
        Function saves the current refractive index map to an hdf5 file:
        ================================================================
        
        Arguments:
        
        - filename: Filename.

        """
        if filename[-3:] == '.h5':
            import h5py
            f = h5py.File(filename,'w')
            map3d = f.create_dataset('data', self._old_map3d.shape, self._old_map3d.dtype)
            map3d[:,:,:] = self._old_map3d[:,:,:]
            f['voxel_dimensions_in_m'] = self._old_map3d_dx
            f.close()
        else:
            log(logger.error,"Invalid filename extension, has to be \'.h5\'.")

    #def load_map3d(self,filename):
    #    """
    #    Function loads refractive index map from an hdf5 file:
    #    ======================================================
    #    
    #    Arguments:
    #    
    #    - filename: Filename.
    #
    #    """
    #
    #    if filename[-3:] == '.h5':
    #        try:
    #            import spimage
    #            img = spimage.sp_image_read(filename,0)
    #            self.map3d = img.image.copy()
    #            spimage.sp_image_free(img)
    #        except:
    #            import h5py
    #            try:
    #                f = h5py.File(filename,'r')
    #                self.map3d = f['data'].value.copy()
    #                #if f['voxel_dimensions_in_m'].value != self.dX: config.OUT.write("WARNING: Sampling of map and setup does not match.")
    #                #self.dX = f['voxel_dimensions_in_m'].value
    #                f.close()
    #            except:
    #                f = h5py.File(filename,'r')
    #                self.map3d = f['real'].value.copy()
    #                f.close()
    #    else:
    #        logger.error("Invalid filename extension, has to be \'.h5\'.")
