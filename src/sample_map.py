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

import sys,numpy,time
import logging
logger = logging.getLogger("Condor")
if "utils" not in sys.path: sys.path.append("utils")
import imgutils,condortools
import utils.nfft
import utils.icosahedron
from scipy import constants
from python_tools.imgtools import array_to_array      

from sample_single import AbstractSampleSingle

class SampleMap(AbstractSampleSingle):

    def __init__(self,**kwargs):
        """
        Function initializes SampleMap object:
        ======================================

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
          ADDITIONAL (icosahedron):
          - geometry_euler_angle_0, geometry_euler_angle_1, geometry_euler_angle_2: Euler angles defining orientation of 5-fold axis in relation to z-axis in grid. [0.0,0.0,0.0]
          ADDITIONAL (spheroid):
          - flattening
          - geometry_euler_angle_0, geometry_euler_angle_1: Euler angles defining orientation of singular axis in relation to z-axis in grid. [0.0,0.0]

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
        AbstractSampleSingle.__init__(self,**kwargs)
        self.map3d_fine = None
        self._map3d = None
        self._dX = None
        self._map3d_fine = None
        self._dX_fine = None

        print "parent",self._parent
        print "kwargs",kwargs
        if "dx_fine" in kwargs:
            self.dX_fine = kwargs["dx_fine"]
        else:
            self.dX_fine = self._parent.detector.get_real_space_resolution_element()/float(kwargs.get("oversampling_fine",1.))/numpy.sqrt(2)

        # Map
        if "geometry" in kwargs:
            if kwargs["geometry"] == "icosahedron":
                if "diameter" not in kwargs:
                    logger.error("Cannot initialize SampleMap instance. diameter is a necessary keyword for geometry=icosahedron.") 
                self.put_icosahedron(self._diameter_mean/2.,**kwargs)
            elif kwargs["geometry"] == "spheroid":
                if "flattening" not in kwargs:
                    logger.error("Cannot initialize SampleMap instance. flattening is a necessary keyword for geometry=spheroid.")
                a = condortools.to_spheroid_semi_diameter_a(self._diameter_mean,kwargs["flattening"])
                c = condortools.to_spheroid_semi_diameter_c(self._diameter_mean,kwargs["flattening"])
                self.put_spheroid(a,c,**kwargs)
            elif kwargs["geometry"] == "sphere":
                self.put_sphere(self._diameter_mean/2.,**kwargs)
            if kwargs["geometry"] == "cube":
                if "edge_length" in kwargs:
                    logger.error("edge_length is a depreciated keyword for geometry=cube. Please replace it by diameter.") 
                self.put_cube(self._diameter_mean,**kwargs)
            elif kwargs["geometry"] == "custom":
                import h5py
                f = h5py.File(kwargs.get("filename","./sample.h5"),"r")
                if len(f.items()) == 1:
                    d = f.items()[0][1][:,:,:]
                else:
                    d = f["data"][:,:,:]
                s = numpy.array(d.shape)
                if not numpy.all(s==s[0]):
                    logger.error("Condor only accepts maps with equal dimensions.")
                    return
                self.map3d_fine = d
                f.close()
                if "dx_fine" in kwargs:
                    self.dX_fine = kwargs["dx_fine"]
                elif "diameter" in kwargs:
                    self.dX_fine = self._diameter_mean/float(s[0])
        else:
            if "map3d_fine" in kwargs:
                s = numpy.array(kwargs["map3d_fine"].shape)
                if not numpy.all(s==s[0]):
                    logger.error("Condor only accepts maps with equal dimensions.")
                    return
                self.map3d_fine = kwargs['map3d_fine']
            else:
                # default
                N = kwargs.get("N_fine",1)
                self.map3d_fine = numpy.zeros(shape=(N,N,N),dtype="float64")

        self._after_init(**kwargs)

    def propagate_single(self,detector0=None,source0=None):
        # scattering amplitude from dn-map: F = F0 DFT{dn} dV
        if source0 == None:
            source = self._parent.source
        else:
            source = source0
        if detector0 == None:
            detector = self._parent.detector
        else:
            detector = detector0

        self.dX = detector.get_real_space_resolution_element() / numpy.sqrt(2)
        map3d = self._get_map3d()
            
        dn_map3d = numpy.array(map3d,dtype="complex128") * self._get_dn()
        self.dn_map3d = dn_map3d

        # scattering vector grid
        q_scaled = detector.generate_qmap(nfft_scaled=True,euler_angle_0=self.euler_angle_0,euler_angle_1=self.euler_angle_1,euler_angle_2=self.euler_angle_2)
        logger.debug("Propagate pattern of %i x %i pixels." % (q_scaled.shape[1],q_scaled.shape[0]))
        q_reshaped = q_scaled.reshape(q_scaled.shape[0]*q_scaled.shape[1],3)
                
        # Check inputs
        invalid_mask = (abs(q_reshaped)>0.5)
        if (invalid_mask).sum() > 0:
            q_reshaped[invalid_mask] = 0.
            logger.debug("%i invalid pixel positions." % invalid_mask.sum())
            
        logger.debug("Map3d input shape: (%i,%i,%i), number of dimensions: %i, sum %f" % (dn_map3d.shape[0],dn_map3d.shape[1],dn_map3d.shape[2],len(list(dn_map3d.shape)),abs(dn_map3d).sum()))
        if (numpy.isfinite(dn_map3d)==False).sum() > 0:
            logger.warning("There are infinite values in the map3d of the object.")
        logger.debug("Scattering vectors shape: (%i,%i); Number of dimensions: %i" % (q_reshaped.shape[0],q_reshaped.shape[1],len(list(q_reshaped.shape))))
        if (numpy.isfinite(q_reshaped)==False).sum() > 0:
            logger.warning("There are infinite values in the scattering vectors.")
        # NFFT
        fourierpattern = utils.nfft.nfft(dn_map3d,q_reshaped)
        # Check output - masking in case of invalid values
        if (invalid_mask).sum() > 0:
            fourierpattern[numpy.any(invalid_mask)] = numpy.nan
        # reshaping
        fourierpattern = numpy.reshape(fourierpattern,(q_scaled.shape[0],q_scaled.shape[1]))

        logger.debug("Got pattern of %i x %i pixels." % (fourierpattern.shape[1],fourierpattern.shape[0]))
        qmap3d = detector.generate_qmap_ori(nfft_scaled=True)

        F = self._get_F0(source,detector) * fourierpattern * self.dX**3
    
        return {"amplitudes": F, "F0": self._get_F0(source, detector),
                "euler_angle_0":self.euler_angle_0,"euler_angle_1":self.euler_angle_1,"euler_angle_2":self.euler_angle_2,
                "dX3": self.dX**3, "grid": q_reshaped, 'qmap3d': qmap3d}
        
    def _get_map3d(self):
        map3d = None
        if self.dX_fine > self.dX:
            logger.error("Finer real space sampling required for chosen geometry.")
            return
        # has map3d_fine the required real space grid?
        if map3d == None and abs(self.dX_fine/self.dX-1) < 0.001:
            # ok, we'll take the fine map
            map3d = self.map3d_fine
            logger.debug("Using the fine map for propagtion.")
            self._map3d = self.map3d_fine
        # do we have an interpolated map?
        if map3d == None and self._dX != None:
            # does it have the right spacing?
            if abs(self._dX/self.dX-1) < 0.001:
                # are the shapes of the original fine map and our current fine map the same?
                if numpy.all(numpy.array(self.map3d_fine.shape)==numpy.array(self._map3d_fine.shape)):
                    # is the grid of the original fine map and the current fine map the same?
                    if self.dX_fine == self._dX_fine:
                        # are the values of the original fine map and the cached fine map the same?
                        if numpy.all(self.map3d_fine==self._map3d_fine):
                            # ok, we take the cached map!
                            map3d = self._map3d
                            logger.debug("Using the cached interpolated map for propagtion.")
        # do we have to do interpolation?
        if map3d == None and self.dX_fine < self.dX:
            from scipy import ndimage
            f = self.dX_fine/self.dX
            N_mapfine = self.map3d_fine.shape[0]
            L_mapfine = (N_mapfine-1)*self.dX_fine
            N_map = int(numpy.floor((N_mapfine-1)*f))+1
            L_map = (N_map-1)*self.dX
            gt = numpy.float64(numpy.indices((N_map,N_map,N_map)))/float(N_map-1)*(N_mapfine-1)*L_map/L_mapfine
            map3d = ndimage.map_coordinates(self.map3d_fine, gt, order=3)
            # Cache interpolated data 
            self._map3d = map3d
            self._dX = N_mapfine/(1.*N_map)*self.dX_fine
            # Cace fine data for later decision whether or not the interpolated map can be used again
            self._map3d_fine = self.map3d_fine
            self._dX_fine = self.dX_fine
            logger.debug("Using a newly interpolated map for propagtion.")
        return map3d


    def put_custom_map(self,map_add,**kwargs):
        unit = kwargs.get("unit","meter")
        p = numpy.array([kwargs.get("z",0.),kwargs.get("y",0.),kwargs.get("x",0.)])
        if unit == "pixel":
            pass
        else:
            p = p/self.dX_fine
        origin = kwargs.get("origin","middle")
        mode = kwargs.get("mode","factor")
        dn = kwargs.get("dn",None)
        if dn == None:
            factor = 1.
        else:
            factor = abs(dn)/abs(self._get_dn())
        if self.map3d_fine == None:
            self.map3d_fine = numpy.array(map_add,dtype="float64")
            return
        else:
            self.map3d_fine = imgtools.array_to_array(map_add,self.map3d_fine,p,origin,mode,0.,factor)

    def put_sphere(self,radius,**kwargs):
        nR = radius/self.dX_fine
        N = int(round((nR*1.2)*2))
        spheremap = make_sphere_map(N,nR)
        self.put_custom_map(spheremap,**kwargs)
 
    def put_spheroid(self,a,c,**kwargs):
        e0 = kwargs.get("geometry_euler_angle_0")
        e1 = kwargs.get("geometry_euler_angle_1")
        e2 = kwargs.get("geometry_euler_angle_2")
        # maximum radius
        Rmax = max([a,c])
        # maximum radius in pixel
        nRmax = Rmax/self.dX_fine
        # dimensions in pixel
        nA = a/self.dX_fine
        nC = c/self.dX_fine
        # leaving a bit of free space around spheroid
        N = int(round((nRmax*1.2)*2))
        spheromap = make_spheroid_map(N,nA,nC,e0,e1,e2)
        self.put_custom_map(spheromap,**kwargs)

    def put_icosahedron(self,radius,**kwargs):
        e0 = kwargs.get("geometry_euler_angle_0")
        e1 = kwargs.get("geometry_euler_angle_1")
        e2 = kwargs.get("geometry_euler_angle_2")
        # icosahedon size parameter
        a = radius*(16*numpy.pi/5.0/(3+numpy.sqrt(5)))**(1/3.0)
        # radius at corners in meter
        Rmax = numpy.sqrt(10.0+2*numpy.sqrt(5))*a/4.0 
        # radius at corners in pixel
        nRmax = Rmax/self.dX_fine 
        # leaving a bit of free space around icosahedron 
        N = int(numpy.ceil(2.3*(nRmax)))
        icomap = make_icosahedron_map(N,nRmax,e0,e1,e2)
        self.put_custom_map(icomap,**kwargs)

    def put_cube(self,a,**kwargs):
        e0 = kwargs.get("geometry_euler_angle_0")
        e1 = kwargs.get("geometry_euler_angle_1")
        e2 = kwargs.get("geometry_euler_angle_2")
        # edge_length in pixels
        nel = a/self.dX_fine 
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
        self.put_custom_map(cubemap,**kwargs)        

    def plot_map3d_fine(self,mode='surface'):
        try:
            from enthought.mayavi import mlab
        except:
            from mayavi import mlab
        if mode=='planes':
            s = mlab.pipeline.scalar_field(abs(self.map3d_fine))
            plane_1 = mlab.pipeline.image_plane_widget(s,plane_orientation='x_axes',
                                                       slice_index=self.map3d_fine.shape[2]/2)
            plane_2 = mlab.pipeline.image_plane_widget(s,plane_orientation='y_axes',
                                                       slice_index=self.map3d_fine.shape[1]/2)
            mlab.show()
        elif mode=='surface':
            mlab.contour3d(abs(self.map3d_fine))
        else:
            logger.error("No valid mode given.")
            
    def plot_fmap3d_fine(self):
        from enthought.mayavi import mlab
        import spimage
        M = spimage.sp_image_alloc(self.map3d_fine.shape[0],self.map3d_fine.shape[1],self.map3d_fine.shape[2])
        M.image[:,:,:] = self.map3d_fine[:,:,:]
        M.mask[:,:,:] = 1
        fM = spimage.sp_image_fftw3(M)
        fM.mask[:,:,:] = 1
        fsM = spimage.sp_image_shift(fM)
        self.fmap3d_fine = abs(fsM.image).copy()
        self.fmap3d_fine[self.fmap3d!=0] = numpy.log10(self.fmap3d_fine[self.fmap3d_fine!=0])
        
        s = mlab.pipeline.scalar_field(self.fmap3d_fine)
        plane_1 = mlab.pipeline.image_plane_widget(s,plane_orientation='x_axes',
                                                   slice_index=self.fmap3d_fine.shape[2]/2)
        plane_2 = mlab.pipeline.image_plane_widget(s,plane_orientation='y_axes',
                                                   slice_index=self.fmap3d_fine.shape[1]/2)
        mlab.show()


    def save_map3d_fine(self,filename):
        """
        Function saves the current refractive index map to an hdf5 file:
        ================================================================
        
        Arguments:
        
        - filename: Filename.

        """
        if filename[-3:] == '.h5':
            import h5py
            f = h5py.File(filename,'w')
            map3d_fine = f.create_dataset('data', self.map3d_fine.shape, self.map3d_fine.dtype)
            map3d_fine[:,:,:] = self.map3d_fine[:,:,:]
            f['voxel_dimensions_in_m'] = self.dX_fine
            f.close()
        else:
            logger.error("Invalid filename extension, has to be \'.h5\'.")

    def load_map3d_fine(self,filename):
        """
        Function loads refractive index map from an hdf5 file:
        ======================================================
        
        Arguments:
        
        - filename: Filename.

        """

        if filename[-3:] == '.h5':
            try:
                import spimage
                img = spimage.sp_image_read(filename,0)
                self.map3d_fine = img.image.copy()
                spimage.sp_image_free(img)
            except:
                import h5py
                try:
                    f = h5py.File(filename,'r')
                    self.map3d_fine = f['data'].value.copy()
                    #if f['voxel_dimensions_in_m'].value != self.dX: config.OUT.write("WARNING: Sampling of map and setup does not match.")
                    #self.dX = f['voxel_dimensions_in_m'].value
                    f.close()
                except:
                    f = h5py.File(filename,'r')
                    self.map3d_fine = f['real'].value.copy()
                    f.close()
        else:
            logger.error("Invalid filename extension, has to be \'.h5\'.")

    def get_area(self,mode='auto'):
        """
        Function returns projected area of the sample along current beam axis:
        ======================================================================

        """
        if self.diameter is not None: return numpy.pi*(self.diameter/2.)**2
        else: return None

def make_icosahedron_map(N,nRmax,euler1=0.,euler2=0.,euler3=0.):
    logger.debug("Building icosahedral geometry")
    logger.debug("Grid: %i x %i x %i (%i voxels)" % (N,N,N,N**3))
    t0 = time.time()
    icomap = utils.icosahedron.icosahedron(N,nRmax,(euler1,euler2,euler3))
    t1 = time.time()
    logger.debug("Built map within %f seconds." % (t1-t0))
    return icomap

def make_icosahedron_map_old(N,nRmax,euler1=0.,euler2=0.,euler3=0.):
    na = nRmax/numpy.sqrt(10.0+2*numpy.sqrt(5))*4.
    nRmin = numpy.sqrt(3)/12*(3.0+numpy.sqrt(5))*na # radius at faces
    logger.debug("Building icosahedral geometry")
    n_list = imgutils.get_icosahedron_normal_vectors(euler1,euler2,euler3)
    X,Y,Z = 1.0*numpy.mgrid[0:N,0:N,0:N]
    X = X - (N-1)/2.
    Y = Y - (N-1)/2.
    Z = Z - (N-1)/2.
    logger.debug("Grid: %i x %i x %i (%i voxels)" % (N,N,N,N**3))
    icomap = numpy.zeros((len(n_list),N,N,N))
    # calculate distance of all voxels to all faces (negative inside, positive outside icosahedron)
    for i in range(len(n_list)):
        icomap[i,:,:,:] = (X*n_list[i][2]+Y*n_list[i][1]+Z*n_list[i][0])+nRmin
    s = 1.
    M = icomap.copy()
    temp = abs(M)<0.5*s
    icomap[temp] = 0.5+icomap[temp]/s
    icomap[M<(-0.5)*s] = 0
    icomap[M>0.5*s] = 1
    icomap = icomap.min(0)
    return icomap

def make_spheroid_map(N,nA,nB,euler0=0.,euler1=0.,euler2=0.):
    X,Y,Z = 1.0*numpy.mgrid[0:N,0:N,0:N]
    X = X-(N-1)/2.
    Y = Y-(N-1)/2.
    Z = Z-(N-1)/2.
    R_sq = X**2+Y**2+Z**2
    e_c = condortools.rotation(numpy.array([0.0,1.0,0.0]),euler0,euler1,euler2)
    d_sq_c = ((Z*e_c[0])+(Y*e_c[1])+(X*e_c[2]))**2
    r_sq_c = abs( R_sq * (1 - (d_sq_c/(R_sq+numpy.finfo("float32").eps))))
    spheroidmap = r_sq_c/nA**2+d_sq_c/nB**2
    spheroidmap[spheroidmap<=1] = 1
    spheroidmap[spheroidmap>1] = 0
    return spheroidmap

def make_sphere_map(N,nR):
    X,Y,Z = 1.0*numpy.mgrid[0:N,0:N,0:N]
    X = X-(N-1)/2.
    Y = Y-(N-1)/2.
    Z = Z-(N-1)/2.
    R = numpy.sqrt(X**2+Y**2+Z**2)
    spheremap = numpy.zeros(shape=R.shape,dtype="float64")
    spheremap[R<=nR] = 1
    spheremap[abs(nR-R)<0.5] = 0.5+0.5*(nR-R[abs(nR-R)<0.5])
    return spheremap
