import sys,numpy,pylab,time,multiprocessing
import logging
logger = logging.getLogger("Propagator")
if "utils" not in sys.path: sys.path.append("utils")
import config,imgutils,proptools

# Pythontools
import gentools,cxitools,imgtools


class Material:
    """
    A class of a sample object.
    Sample material.

    """
    def __init__(self,parent=None,**kwargs):
        self._parent = parent
        if "massdensity" in kwargs:
            self.materialtype = "custom"
            for key in kwargs:
                if key[0] == 'c' or key == 'massdensity':
                    exec "self." + key + " = args[key]"
                else:
                    logger.error("%s is no valid argument for custom initialization of Material." % key)
                    return
        elif "material_type" in kwargs:
            self.material_type = kwargs['material_type']
            self.massdensity = config.DICT_massdensity[self.material_type]
            self.cH = config.DICT_atomic_composition[self.material_type][0]
            self.cC = config.DICT_atomic_composition[self.material_type][1]
            self.cN = config.DICT_atomic_composition[self.material_type][2]
            self.cO = config.DICT_atomic_composition[self.material_type][3]
            self.cP = config.DICT_atomic_composition[self.material_type][4]
            self.cS = config.DICT_atomic_composition[self.material_type][5]
        else:
            logger.error("No valid arguments for Material initialization.")
            return

    def get_fX(self,element,photon_energy_eV=None):
        """
        Get the scattering factor for an element through linear interpolation.
        """
        if not photon_energy_eV:
            photon_energy_eV = self._parent._parent.source.photon.get_energy("eV")
        SF_X = config.DICT_scattering_factors[element]
        e = config.DICT_physical_constants['e']
        c = config.DICT_physical_constants['c']
        h = config.DICT_physical_constants['h']
        f1 = numpy.interp(photon_energy_eV,SF_X[:,0],SF_X[:,1])
        f2 = numpy.interp(photon_energy_eV,SF_X[:,0],SF_X[:,2])
        return complex(f1,f2) 
 
    def get_n(self,photon_energy_eV=None):
        """
        Obtains complex refractive index.
        Henke (1994): n = 1 - r_0/(2pi) lambda^2 sum_q rho_q f_q(0)
        r_0: classical electron radius
        rho_q: atomic number density of atom species q
        f_q(0): atomic scattering factor (forward scattering) of atom species q
        """

        re = config.DICT_physical_constants['re']
        h = config.DICT_physical_constants['h']
        c = config.DICT_physical_constants['c']
        qe = config.DICT_physical_constants['e']

        if not photon_energy_eV:
            photon_energy_eV = self._parent._parent.source.photon.get_energy("eV")
        photon_wavelength = h*c/photon_energy_eV/qe

        f = self.get_f(photon_energy_eV)
        atom_density = self.get_atom_density()
        
        n = 1 - re/2/numpy.pi * photon_wavelength**2 * f * atom_density

        return n

    def get_dn(self,photon_energy=None):
        return (1-self.get_n(photon_energy))

    # convenience functions
    # n = 1 - delta - i beta
    def get_delta(self,photon_energy_eV=None):
        return (1-self.get_n(photon_energy_eV).real)
    def get_beta(self,photon_energy_eV=None):
        return (-self.get_n(photon_energy_eV).imag)
    def get_photoabsorption_cross_section(self,photon_energy_eV=None):
        re = config.DICT_physical_constants['re']
        h = config.DICT_physical_constants['h']
        c = config.DICT_physical_constants['c']
        qe = config.DICT_physical_constants['e']
        if not photon_energy_eV:
            photon_energy_eV = self._parent._parent.source.photon.get_energy("eV")
        photon_wavelength = h*c/photon_energy_eV/qe
        mu = 2*re*photon_wavelength*self.get_f(photon_energy_eV).imag
        return mu
    def get_transmission(self,thickness,photon_energy_eV=None):
        n = self.get_n(photon_energy_eV)
        mu = self.get_photoabsorption_cross_section(photon_energy_eV)
        rho = self.get_atom_density()
        return numpy.exp(-rho*mu*thickness)

    def get_f(self,photon_energy_eV=None):

        h = config.DICT_physical_constants['h']
        c = config.DICT_physical_constants['c']
        qe = config.DICT_physical_constants['e']

        if not photon_energy_eV:
            photon_energy_eV = self._parent._parent.source.photon.get_energy("eV")
        photon_wavelength = h*c/photon_energy_eV/qe

        atomic_composition = self.get_atomic_composition_dict()

        f_sum = 0
        for element in atomic_composition.keys():
            # sum up average atom factor
            f = self.get_fX(element,photon_energy_eV)
            f_sum += atomic_composition[element] * f
        
        return f_sum


    def get_atom_density(self):
                
        u = config.DICT_physical_constants['u']

        atomic_composition = self.get_atomic_composition_dict()

        M = 0
        for element in atomic_composition.keys():
            # sum up mass
            M += atomic_composition[element]*config.DICT_atomic_mass[element]*u

        number_density = self.massdensity/M
        
        return number_density


    def get_electron_density(self):

        u = config.DICT_physical_constants['u']

        atomic_composition = self.get_atomic_composition_dict()

        M = 0
        Q = 0
        for element in atomic_composition.keys():
            # sum up electrons
            M += atomic_composition[element]*config.DICT_atomic_mass[element]*u
            Q += atomic_composition[element]*config.DICT_atomic_number[element]

        electron_density = Q*self.massdensity/M
        
        return electron_density
        
        
    def get_atomic_composition_dict(self):

        atomic_composition = {}
        
        for key in self.__dict__.keys():
            if key[0] == 'c':
                exec "c_tmp = self." + key
                atomic_composition[key[1:]] = c_tmp 
 
        tmp_sum = float(sum(atomic_composition.values()))
        for element in atomic_composition.keys():
            atomic_composition[element] /= tmp_sum 
        
        return atomic_composition

class Sample:
    def __init__(self,**kwargs):
        self._parent = kwargs.get('parent',None)
        # Material
        materialargs = {}
        if 'massdensity' in kwargs:
            materialargs['massdensity'] = kwargs['massdensity']
            for key in kwargs.keys():
                if key[0] == 'c': materialargs[key] = kwargs[key]
            self.material = Material(self,**materialargs)
        elif "material_type" in kwargs:
            materialargs['material_type'] = kwargs['material_type']
            self.material = Material(self,**materialargs)
        else:
            self.material = None

    def set_random_orientation(self):
        [e0,e1,e2] = proptools.random_euler_angles()
        self.euler_angle_0 = e0
        self.euler_angle_1 = e1
        self.euler_angle_2 = e2

class SampleSphere(Sample):
    """
    A class of the input-object.
    Sample is a homogeneous sphere defined by a radius and a material object.

    """

    def __init__(self,**kwargs):
        Sample.__init__(self,**kwargs)
        reqk = ["size"]
        for k in reqk:
            if k not in kwargs.keys():
                logger.error("Cannot initialize SampleSphere instance. %s is a necessary keyword." % k)
                return
        self.radius = kwargs['size']/2.
        self._parent = kwargs.get('parent',None)

        material_kwargs = kwargs.copy()
        non_material_arguments = ['parent','size']
        for non_material_argument in non_material_arguments:
            try: material_kwargs.pop(non_material_argument)
            except: pass
        material_kwargs['parent'] = self
        self.material = Material(**material_kwargs)

    def get_area(self):
        """ Calculates area of projected sphere """
        return numpy.pi*self.radius**2

class SampleSpheroid(Sample):
    """
    A class of the input-object.
    Sample is a homogeneous spheroid defined by the orthogonal diameters a and c (c is the one along the rotation axis) and a material object.

    """

    def __init__(self,**kwargs):
        Sample.__init__(self,**kwargs)
        reqk = ["diameter_a","diameter_c","theta","phi"]
        for k in reqk:
            if k not in kwargs.keys():
                logger.error("Cannot initialize SampleSpheroid instance. %s is a necessary keyword." % k)
                return
        self.a = kwargs['diameter_a']/2.
        self.c = kwargs['diameter_c']/2.
        self.theta = kwargs["theta"]
        self.phi = kwargs["phi"]

        material_kwargs = kwargs.copy()
        non_material_arguments = ['parent','diameter_a','diameter_c','theta','phi']
        for non_material_argument in non_material_arguments:
            try: material_kwargs.pop(non_material_argument)
            except: pass
        material_kwargs['parent'] = self
        self.material = Material(**material_kwargs)

    def get_area(self):
        """
        Calculates area of projected spheroid
        """
        logger.warning("Calculates area of WRONGLY projected spheroid, fix when there is time.")
        return (4/3.*numpy.pi*self.a**2*self.c)**(2/3.)


class SampleMap(Sample):
    """
    A class of the input-object.
    Class of the sample described by a refractive index map (self.map3d_fine: dn = n-1) 

    """

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
          - dX_fine: Real space sampling distance.
          OR:
          - oversampling_fine: Real space oversampling in relation to detector/source configuration. [1]
        
        OR:
        - geometry: Geometry that shall be generated (icosahedron, sphere, spheroid)
          - oversampling_fine: Additional oversampling for the initial map self.map3d_fine [1.]
          ADDITIONAL (icosahedron):
          - diameter in meter (sphere-volume equivalent diameter)
          - geometry_euler_angle_0, geometry_euler_angle_1, geometry_euler_angle_2: Euler angles defining orientation of 5-fold axis in relation to z-axis in grid. [0.0,0.0,0.0]
          - oversampling_fine: Real space oversampling in relation to detector/source configuration.
          ADDITIONAL (sphere):
          - diameter in meter
          - oversampling_fine: Real space oversampling in relation to detector/source configuration.
          ADDITIONAL (spheroid):
          - diameter_c: diameter along the singular axis (rotation axis of ellipsoid:
          - diameter_a: diameter along orthogonal axis to singular axis
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

        Sample.__init__(self,**kwargs)
        self.euler_angle_0 = kwargs.get("euler_angle_0",0.)
        self.euler_angle_1 = kwargs.get("euler_angle_1",0.)
        self.euler_angle_2 = kwargs.get("euler_angle_2",0.)
        self.map3d_fine = None
        self._map3d = None
        self._dX = None
        self.radius = kwargs.get('radius',None)

        if "dX_fine" in kwargs:
            self.dX_fine = kwargs["dX_fine"]
        elif "oversampling_fine":
            self.dX_fine = self._parent.get_real_space_resolution_element()/float(kwargs["oversampling_fine"])

        # Map
        if "geometry" in kwargs:
            if "geometry" in kwargs:
                if kwargs["geometry"] == "icosahedron":
                    if "diameter" not in kwargs:
                        logger.error("Cannot initialize SampleMap instance. diameter is a necessary keyword for geometry=icosahedron.") 
                    self.put_icosahedron(kwargs["diameter"]/2.,**kwargs)
                    self.radius = kwargs["diameter"]/2.
                elif kwargs["geometry"] == "spheroid":
                    if "diameter_a" not in kwargs or "diameter_c" in kwargs:
                        logger.error("Cannot initialize SampleMap instance. a_diameter and c_diameter are necessary keywords for geometry=spheroid.")
                    self.put_spheroid(kwargs["diameter_a"],kwargs["diameter_c"],k)
                    self.radius = (2*kwargs["diameter_a"]+kwargs["diameter_c"])/3./2.
                elif kwargs["geometry"] == "sphere":
                    if "diameter" not in kwargs:
                        logger.error("Cannot initialize SampleMap instance. diameter is a necessary keyword for geometry=sphere.")                
                    self.put_sphere(self.radius,**kwargs)
                    self.radius = kwargs["diameter"]/2.
        else:
            if "map3d_fine" in kwargs:
                s = numpy.array(self.maprgs["map3d_fine"].shape)
                if not numpy.all(s==s[0]):
                    logger.error("Propagator only accepts maps with equal dimensions.")
                    return
                self.map3d_fine = kwargs['map3d_fine']
            else:
                # default
                N = kwargs.get("N_fine",1)
                self.map3d_fine = numpy.zeros(shape=(N,N,N),dtype="complex64")

    def get_refractive_index_map3d_fine(self):
        if self.material == None:
            dn = 1.
        else:
            dn = self.material.get_dn()
        return self.map3d_fine*dn
        
    def put_custom_map(self,map_add,**kwargs):
        unit = kwargs.get("unit","meter")
        p = numpy.array([kwargs.get("z",0.),kwargs.get("y",0.),kwargs.get("x",0.)])
        if unit == "pixel":
            pass
        else:
            p = p/self.dX_fine
        origin = kwargs.get("origin","middle")
        mode = kwargs.get("mode","sum")
        if self.map3d_fine == None:
            self.map3d_fine = numpy.array(map_add,dtype="complex64")
            return
        else:
            p = numpy.array([z/self.dX_fine,y/self.dX_fine,z/self.dX_fine])
            self.map3d_fine = imgtools.array_to_array(self.map3d_fine,map_add,p,origin,mode)

    def put_sphere(self,radius,**kwargs):
        R_N = radius/self.dX_fine
        size = int(round((R_N*1.2)*2))
        X,Y,Z = 1.0*numpy.mgrid[0:size,0:size,0:size]
        for J in [X,Y,Z]: J = J - size/2.0-0.5
        R = numpy.sqrt(X**2+Y**2+Z**2)
        spheremap = numpy.zeros(shape=R.shape,dtype="complex64")
        spheremap[R<R_N] = 1
        spheremap[abs(R_N-R)<0.5] = 0.5+0.5*(R_N-R[abs(R_N-R)<0.5])
        self.put_custom_map(spheremap,**kwargs)
 
    def put_spheroid(self,a,b,**kwargs):
        e0 = kwargs.get("geometry_euler_angle_0")
        e1 = kwargs.get("geometry_euler_angle_1")
        e2 = kwargs.get("geometry_euler_angle_2")
        # mxaximum radius
        Rmax = max([a,b])
        # maximum radius in pixel
        nRmax = Rmax/self.dX_fine
        # dimensions in pixel
        nA = a/self.dX_fine
        nB = b/self.dX_fine
        # leaving a bit of free space around spheroid
        N = int(round((R_N*1.2)*2))
        spheromap = make_spheroid_map(N,nA,nB,e0,e1,e2)
        self.put_custom_map(spheroidmap,**kwargs)

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

    def project(self,eul_ang0=None,eul_ang1=None,eul_ang2=None,**kwargs):
        """ Projection of 3-dimensional map."""
        if not eul_ang0: eul_ang0 = self.grid_euler_angle_0
        if not eul_ang1: eul_ang1 = self.grid_euler_angle_1
        if not eul_ang2: eul_ang2 = self.grid_euler_angle_2

        if eul_ang0 == 0.0 and eul_ang1 == 0.0 and eul_ang2 == 0.0:
            map2d = numpy.zeros((self.map3d_fine.shape[1],self.map3d_fine.shape[2]))
            for iy in numpy.arange(0,map2d.shape[0]):
                for ix in numpy.arange(0,map2d.shape[1]):
                    map2d[iy,ix] = self.map3d_fine[:,iy,ix].real.sum()*self.dX_fine*numpy.exp(-self.map3d_fine[:,iy,ix].imag.sum()*self.dX_fine)
        else:
            x_0 = proptools.rotation(numpy.array([0.0,0.0,1.0]),eul_ang0,eul_ang1,eul_ang2)
            y_0 = proptools.rotation(numpy.array([0.0,1.0,0.0]),eul_ang0,eul_ang1,eul_ang2)
            z_0 = proptools.rotation(numpy.array([1.0,0.0,0.0]),eul_ang0,eul_ang1,eul_ang2)
            N = self.map3d_fine.shape[0]
            extra_space = 10
            map2d = numpy.zeros(shape=(N,N),dtype="complex64")
            def valid_coordinate(v):
                for i in range(0,3):
                    if v[i] < N and v[i] > 0:
                        pass
                    else:
                        return False
                return True
            for y in range(0,N):
                for x in range(0,N):
                    for z in range(0,N):
                        vector = numpy.array([(N-1)/2.0,(N-1)/2.0,(N-1)/2.0])-(N-1)/2.0*(x_0+y_0+z_0)+(z*z_0+y*y_0+x*x_0)
                        cases = [numpy.floor,lambda x: 1.0 + numpy.floor(x)]
                        for y_func in cases:
                            for x_func in cases:
                                for z_func in cases:
                                    if valid_coordinate(numpy.array([z_func(vector[0]),y_func(vector[1]),x_func(vector[2])])):
                                        map2d[y,x] +=\
                                            self.dX_fine*\
                                            (1.0-abs(z_func(vector[0]) - vector[0]))*\
                                            (1.0-abs(y_func(vector[1]) - vector[1]))*\
                                            (1.0-abs(x_func(vector[2]) - vector[2]))*\
                                            self.map3d_fine[int(z_func(vector[0])),
                                                  int(y_func(vector[1])),
                                                  int(x_func(vector[2]))]
        if "padding" in kwargs.keys():
            if kwargs["padding"]=="squared":
                N_new = max([map2d.shape[0],map2d.shape[1]])
                map2d_new = numpy.zeros(shape=(N_new,N_new))
                map2d_new[(N_new-map2d.shape[0])/2:(N_new-map2d.shape[0])/2+map2d.shape[0],
                          (N_new-map2d.shape[1])/2:(N_new-map2d.shape[1])/2+map2d.shape[1]] = map2d[:,:]
                map2d = map2d_new
        return map2d

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
        if mode == 'auto':
            if self.radius != None: return numpy.pi*self.radius**2
        projection = self.project(self.grid_euler_angle_0,self.grid_euler_angle_1,self.grid_euler_angle_2)
        binary = numpy.ones(shape=projection.shape)
        binary[abs(projection)<numpy.median(abs(projection))] = 0
        area = binary.sum()*self.dX_fine**2
        return area

    
    

    
def make_icosahedron_map(N,nRmax,euler1=0.,euler2=0.,euler3=0.,s=1.):

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

    M = icomap.copy()
    temp = abs(M)<0.5*s
    icomap[temp] = 0.5+icomap[temp]/s
    icomap[M<(-0.5)*s] = 0
    icomap[M>0.5*s] = 1
    icomap = icomap.min(0)

    return icomap


def make_spheroid_map(N,nA,nB,euler0=0.,euler1=0.,euler2=0.,s=1.):
    X,Y,Z = 1.0*numpy.mgrid[0:N,0:N,0:N]
    for J in [X,Y,Z]: J -= N/2.0-0.5        
    R_sq = X**2+Y**2+Z**2
    e_c = proptools.rotation(numpy.array([0.0,0.0,1.0]),e0,e1,e2)
    d_sq_c = ((Z*e_c[0])+(Y*e_c[1])+(X*e_c[2]))**2
    r_sq_c = abs( R_sq * (1 - (d_sq_c/R_sq)))
    spheroidmap = r_sq_c/a_N**2+d_sq_c/b_N**2
    spheroidmap[spheroidmap<=1] = 1
    spheroidmap[spheroidmap>1] = 0
    logger.info("Smoothing by a factor of %f\n" % s)
    spheroidmap = abs(imgutils.smooth3d(spheroidmap,s))
    spheroidmap /= spheroidmap.max()
    return spheroidmap
