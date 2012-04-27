import sys,pylab,time,multiprocessing
if "utils" not in sys.path: sys.path.append("utils")
import config,imgutils,tools

class Material:
    """
    A class of a sample object.
    Sample material.

    """
    def __init__(self,parent=None,**args):
        self._parent = parent
        keys = args.keys()
        try:
            args['massdensity']
            for key in keys:
                if key[0] == 'c' or key == 'massdensity':
                    exec "self." + key + " = args[key]"
                else:
                    print "ERROR: %s is no valid argument for initialization of Material." % key
        except:
            if keys == []:
                materialtype = 'protein'
            else:
                materialtype = args['materialtype']
            self.massdensity = config.DICT_massdensity[materialtype]
            self.cH = config.DICT_atomic_composition[materialtype][0]
            self.cC = config.DICT_atomic_composition[materialtype][1]
            self.cN = config.DICT_atomic_composition[materialtype][2]
            self.cO = config.DICT_atomic_composition[materialtype][3]
            self.cP = config.DICT_atomic_composition[materialtype][4]
            self.cS = config.DICT_atomic_composition[materialtype][5]
    
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
        f1 = pylab.interp(photon_energy_eV,SF_X[:,0],SF_X[:,1])
        f2 = pylab.interp(photon_energy_eV,SF_X[:,0],SF_X[:,2])
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
            
        n = 1 - re/2/pylab.pi * photon_wavelength**2 * f

        return n


    def get_f(self,photon_energy_eV=None):

        h = config.DICT_physical_constants['h']
        c = config.DICT_physical_constants['c']
        qe = config.DICT_physical_constants['e']
        
        atom_density = self.get_atom_density()

        if not photon_energy_eV:
            photon_energy_eV = self._parent._parent.source.photon.get_energy("eV")
        photon_wavelength = h*c/photon_energy_eV/qe

        atomic_composition = self.get_atomic_composition_dict()

        rhof_sum = 0
        for element in atomic_composition.keys():
            # sum up average atom factor
            f = self.get_fX(element,photon_energy_eV)
            rhof_sum += atomic_composition[element] * atom_density * f

        return rhof_sum


    def get_atom_density(self):
                
        u = config.DICT_physical_constants['u']

        atomic_composition = self.get_atomic_composition_dict()

        M = 0
        for element in atomic_composition.keys():
            # sum up average atom density
            M += atomic_composition[element]*config.DICT_atomic_mass[element]*u

        number_density = self.massdensity/M
        
        return number_density

        
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



class SampleMap:
    """
    A class of the input-object.
    Class of the sample described by a refractive index map (map3d: dn = n-1) 

    """

    def __init__(self,**kwargs):
        """
        Function initializes SampleMap object:
        ======================================

        Arguments:

        Keyword arguments (if not given variable \'X\' is set to default value \'[X_default]\'):
        
        EITHER (default):
        - N: Edge length in number of voxels [1]
        - dX Real space sampling distance. [Taken from input object configuration]

        OR:
        - map3d: Cubic 3d refractive index map.
        
        ADDITIONAL:
        - euler_angle_0, euler_angle_1, euler_angle_2: Euler angles defining orientation of particle in the beam. [0.0,0.0,0.0]
        - radius: Characteristic radius of the object. [None]
        - parent: Input object to which this SampleMap object shall be linked. [None]

        """
        if 'map3d' in kwargs: 
            if kwargs['map3d'].shape[0] != kwargs['map3d'].shape[1] \
                    or kwargs['map3d'].shape[0] != kwargs['map3d'].shape[2] \
                    or kwargs['map3d'].shape[1] != kwargs['map3d'].shape[2]:
                print "ERROR: Only accept cubic maps."
                return
            self.map3d = kwargs['map3d']
        else:
            N = kwargs.get('N',1)
            self.map3d = pylab.zeros(shape=(N,N,N))
        self._parent = kwargs.get('parent',None)
        self.radius = kwargs.get('radius',None)
        self.dX = kwargs.get('dX',self._parent.get_real_space_resolution_element()/(1.0*self._parent.propagation.rs_oversampling))
        self.euler_angle_1 = kwargs.get('euler_angle_1',0.0)
        self.euler_angle_2 = kwargs.get('euler_angle_2',0.0)
        self.euler_angle_3 = kwargs.get('euler_angle_3',0.0)

    def set_random_orientation(self):
        [e1,e2,e3] = tools.random_euler_angles()
        self.euler_angle_1 = e1
        self.euler_angle_2 = e2
        self.euler_angle_3 = e3
        
    def put_custom_map(self,map_add,x,y,z):
        origin = pylab.array([(self.map3d.shape[2]-1)/2.0,
                              (self.map3d.shape[1]-1)/2.0,
                              (self.map3d.shape[0]-1)/2.0])
        zmin = round(origin[2]+z-map_add.shape[0]/2.0)
        ymin = round(origin[1]+y-map_add.shape[1]/2.0)
        xmin = round(origin[0]+x-map_add.shape[2]/2.0)
        self.map3d[zmin:zmin+map_add.shape[0],
                   ymin:ymin+map_add.shape[1],
                   xmin:xmin+map_add.shape[2]] += map_add[:,:,:]

    # TESTING REQUIRED
    def smooth_map3d(self,factor):
        self.map3d = imgutils.smooth3d(self.map3d,factor)

    def put_sphere(self,radius,x=0,y=0,z=0,**materialargs):
        """
        Function adds densitymap of homogeneous sphere to 3-dimensional densitymap:
        ===========================================================================

        Arguments:
       
        - radius: Radius of the sphere in m. [no default value]
        - x,y,z: center position (measured from current origin (0,0,0)) in m. [0.,0.,0.]

        Keyword arguments:

        EITHER:
        - materialtype: Predefined materialtype ('protein','virus','cell','latexball','water',... For a complete list look up \'config.py\'.)
        
        OR:
        - c\'X\': fraction of element \'X\' of atomic number density in material (normalization is done internally).
        - massdensity: Massdensity of material in kg/m^3

        """
        material_obj = Material(self,**materialargs)
        dn = 1.0-material_obj.get_n()
        R_N = radius/self.dX
        size = int(round((R_N*1.2)*2))
        X,Y,Z = 1.0*pylab.mgrid[0:size,0:size,0:size]
        for J in [X,Y,Z]: J -= size/2.0-0.5
        R = pylab.sqrt(X**2+Y**2+Z**2)
        spheremap = pylab.zeros_like(R)
        spheremap[R<R_N] = 1
        spheremap[abs(R_N-R)<0.5] = 0.5+0.5*(R_N-R[abs(R_N-R)<0.5])
        spheremap *= dn
        if self.map3d.max() == 0 and not x and not y and not z:
            self.map3d = spheremap
        else:
            x_N = x/self.dX
            y_N = y/self.dX
            z_N = z/self.dX
            self.put_custom_map(spheremap,x_N,y_N,z_N)
 
    def put_spheroid(self,a,b,x=None,y=None,z=None,eul_ang1=0.0,eul_ang2=0.0,eul_ang3=0.0,**materialargs):
        """
        Function add densitymap of homogeneous ellipsoid to 3-dimensional densitymap:
        =============================================================================

        Arguments:
        
        - a,b: Radii in meter.
        - x,y,z: Center coordinats (from current origin) [middle]

        Keyword arguments:

        EITHER:
        - materialtype: Predefined materialtype ('protein','virus','cell','latexball','water',... For a complete list look up \'config.py\'.)
        
        OR:
        - c\'X\': fraction of element \'X\' of atomic number density in material (normalization is done internally).
        - massdensity: Massdensity of material in kg/m^3

        """

        material_obj = Material(self,**materialargs)
        dn = 1.0-material_obj.get_n()        
        R_N = max([a,b])/self.dX
        a_N = a/self.dX
        b_N = b/self.dX
        size = int(round((R_N*1.2)*2))
        X,Y,Z = 1.0*pylab.mgrid[0:size,0:size,0:size]
        for J in [X,Y,Z]: J -= size/2.0-0.5
        
        if eul_ang1 != 0.0 or eul_ang2 != 0.0 or eul_ang3 != 0.0:
            def rotate_X(v,alpha):
                rotM = pylab.array([[1,0,0],[0,pylab.cos(alpha),-pylab.sin(alpha)],[0,pylab.sin(alpha),pylab.cos(alpha)]])
                return pylab.dot(rotM,v)
            def rotate_Z(v,alpha):
                rotM = pylab.array([[pylab.cos(alpha),-pylab.sin(alpha),0],[pylab.sin(alpha),pylab.cos(alpha),0],[0,0,1]])
                return pylab.dot(rotM,v)
            def do_basic_rotation(v,e1,e2,e3):
                new = rotate_Z(v,e1)
                new = rotate_X(new,e2)
                new = rotate_Z(new,e3)
                return new
            print do_basic_rotation(pylab.array([2.3,2.3,2.3]),0.1,0.1,0.1)
            for xi in pylab.arange(0,size,1.0):
                for yi in pylab.arange(0,size,1.0):
                    for zi in pylab.arange(0,size,1.0):
                        new = pylab.array([Z[zi,yi,xi],Y[zi,yi,xi],X[zi,yi,xi]])
                        new = do_basic_rotation(new,eul_ang1,eul_ang2,eul_ang3)
                        X[zi,yi,xi] = new[2]
                        Y[zi,yi,xi] = new[1]
                        Z[zi,yi,xi] = new[0]
        spheroidmap = (X**2+Y**2)/a_N**2+Z**2/b_N**2
        spheroidmap[spheroidmap<=1] = 1
        spheroidmap[spheroidmap>1] = 0
        spheroidmap *= dn
        if self.map3d.max() == 0 and not x and not y and not z:
            self.map3d = spheroidmap
        else:
            if x == None and y == None and z == None:
                x=round(self.map3d.shape[2]/2.0)-0.5
                y=round(self.map3d.shape[1]/2.0)-0.5
                z=round(self.map3d.shape[0]/2.0)-0.5
            x_N = x/self.dX
            y_N = y/self.dX
            z_N = z/self.dX
            self.put_custom_map(spheroidmap,x_N,y_N,z_N)
   
    def put_goldball(self,radius,x=None,y=None,z=None):
        self.put_sphere(radius,x,y,z,cAu=1,massdensity=config.DICT_massdensity['Au'])        

    def put_icosahedral_virus(self,radius,eul_ang1=0.0,eul_ang2=0.0,eul_ang3=0.0,x=0.,y=0.,z=0.):
        dn = self._makedm_icosahedron(radius,eul_ang1,eul_ang2,eul_ang3,materialtype="virus")
        if self.map3d.max() == 0.0 and x==0. and y==0. and z==0.:
            self.map3d = dn
        else:
            x_N = x/self.dX
            y_N = y/self.dX
            z_N = z/self.dX
            self.put_custom_map(dn,x_N,y_N,z_N)


    def _makedm_icosahedron(self,radius,eul_ang1=0.0,eul_ang2=0.0,eul_ang3=0.0,**materialargs):  
        """
        Function returns a refractive index map of a homogeneous icosahedron:
        =====================================================================

        Arguments:

        - radius: Radius in meter (volume equals volume of a sphere with given radius). [no default value]
        - eul_ang1,eul_ang2,eul_ang3: orientation defined by Euler angles in rad. [0.0,0.0,0.0]
        
        Keyword arguments:

        EITHER:
        - materialtype: Predefined materialtype ('protein','virus','cell','latexball','water',... For a complete list look up \'config.py\'.)
        
        OR:
        - c\'X\': fraction of element \'X\' of atomic number density in material (normalization is done internally).
        - massdensity: Massdensity of material in kg/m^3

        """

        t_start = time.time()

        a = radius*(16*pylab.pi/5.0/(3+pylab.sqrt(5)))**(1/3.0)
        Rmax = pylab.sqrt(10.0+2*pylab.sqrt(5))*a/4.0 # radius at corners
        Rmin = pylab.sqrt(3)/12*(3.0+pylab.sqrt(5))*a # radius at faces
        nRmax = Rmax/self.dX
        nRmin = Rmin/self.dX
        # smooth
        s = 1.0
        N = int(pylab.ceil(2*(nRmax+pylab.ceil(s))))
        r_pix = self.dX*(3/(4*pylab.pi))**(1/3.0)
        
        config.OUT.write("... build icosahedral geometry ...\n")
        
        # construct normal vectors of faces
        phi = (1+pylab.sqrt(5))/2.0
        # normal vectors for every vertice
        x1 = pylab.array([0.0,1.0,phi])
        x2 = pylab.array([0.0,1.0,-phi])
        x3 = pylab.array([0.0,-1.0,phi])
        x4 = pylab.array([0.0,-1.0,-phi]) 
        x5 = pylab.array([1.0,phi,0.0])
        x6 = pylab.array([1.0,-phi,0.0])
        x7 = pylab.array([-1.0,phi,0.0])
        x8 = pylab.array([-1.0,-phi,0.0])
        x9 = pylab.array([phi,0.0,1.0])
        x10 = pylab.array([-phi,0.0,1.0])
        x11 = pylab.array([phi,0.0,-1.0])
        x12 = pylab.array([-phi,0.0,-1.0])
        X = [x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12]
        def cont_element(el,l):
            for i in range(0,len(l)):
                if (el == l[i]).all():
                    return True
            return False
        an = round(pylab.dot(x5,x1))
        def angles_match(y1,y2,y3):
            if round(pylab.dot(y1,y2)) == an and round(pylab.dot(y2,y3)) == an and round(pylab.dot(y3,y1)) == an:
                return True
            else:
                return False
        n_list = []
        for i in range(0,len(X)):
            for j in range(0,len(X)):
                for k in range(0,len(X)):
                    n = (X[i]+X[j]+X[k])/6*a/self.dX
                    if angles_match(X[i],X[j],X[k]) and not cont_element(n,n_list):
                        n_list.append(n)

        config.OUT.write("... rotate coordinate system in respect to given Euler angles ...\n")        
        def rotate_X(v,alpha):
            rotM = pylab.array([[1,0,0],[0,pylab.cos(alpha),-pylab.sin(alpha)],[0,pylab.sin(alpha),pylab.cos(alpha)]])
            return pylab.dot(rotM,v)
        def rotate_Z(v,alpha):
            rotM = pylab.array([[pylab.cos(alpha),-pylab.sin(alpha),0],[pylab.sin(alpha),pylab.cos(alpha),0],[0,0,1]])
            return pylab.dot(rotM,v)
        for i in range(0,len(n_list)):
            n_list[i] = rotate_Z(n_list[i],eul_ang1)
            n_list[i] = rotate_X(n_list[i],eul_ang2)
            n_list[i] = rotate_Z(n_list[i],eul_ang3)

        #t_1 = time.time()        
        config.OUT.write("... %i x %i x %i grid (%i voxels) ...\n" % (N,N,N,N**3))
        icomap = pylab.ones((N,N,N))
        positions = []
        for iz in range(0,N):
            for iy in range(0,N):
                for ix in range(0,N):
                    positions.append([iz,iy,ix])
        positions = pylab.array(positions)

        #t_2 = time.time()
        pool = multiprocessing.Pool()
        positions_pool = []
        results_pool = []
        N_chunks = multiprocessing.cpu_count()#int(pylab.ceil(len(positions)/(1.0*chunksize)))
        chunksize = int(pylab.ceil(len(positions)/(1.0*N_chunks)))

        for chunk in range(0,N_chunks):
            config.OUT.write("Starting process %i / %i\n" % ((chunk+1),N_chunks))
            if chunk < (N_chunks-1): positions_pool.append(positions[chunksize*chunk:chunksize*(chunk+1),:])
            else: positions_pool.append(positions[chunksize*chunk:,:])
            results_pool.append(pool.apply_async(imgutils.cut_edges,(positions_pool[-1],n_list,radius,self.dX)))
        pool.close()
        pool.join()
        cuts = pylab.ones(len(positions))
        for chunk in range(0,N_chunks):
            config.OUT.write("Receiving results from process %i / %i\n" % ((chunk+1),N_chunks))
            res = results_pool[chunk].get()
            for i in range(0,len(positions_pool[chunk])):
                [iz,iy,ix] = positions_pool[chunk][i]
                icomap[iz,iy,ix] = res[i]
        t_stop = time.time()
        config.OUT.write("Time: %f sec\n" % (t_stop-t_start))

        material_obj = Material(self,**materialargs)
        dn = 1.0 - material_obj.get_n()

        return dn*icomap
    

    def project(self,eul_ang1=None,eul_ang2=None,eul_ang3=None,**kwargs):
        """ Projection of 3-dimensional map."""
        if not eul_ang1: eul_ang1 = self.euler_angle_1
        if not eul_ang2: eul_ang2 = self.euler_angle_2
        if not eul_ang3: eul_ang3 = self.euler_angle_3

        if eul_ang1*eul_ang2*eul_ang3 == 0.0:
            map2d = pylab.zeros((self.map3d.shape[1],self.map3d.shape[2]))
            for iy in pylab.arange(0,map2d.shape[0]):
                for ix in pylab.arange(0,map2d.shape[1]):
                    map2d[iy,ix] = self.map3d[:,iy,ix].real.sum()*self.dX*pylab.exp(-self.map3d[:,iy,ix].imag.sum()*self.dX)
        else:
            def rotate_X(v,alpha):
                rotM = pylab.array([[1,0,0],[0,pylab.cos(alpha),-pylab.sin(alpha)],[0,pylab.sin(alpha),pylab.cos(alpha)]])
                return pylab.dot(rotM,v)
            def rotate_Z(v,alpha):
                rotM = pylab.array([[pylab.cos(alpha),-pylab.sin(alpha),0],[pylab.sin(alpha),pylab.cos(alpha),0],[0,0,1]])
                return pylab.dot(rotM,v)
            x_0 = pylab.array([0.0,0.0,1.0])
            x_0 = rotate_Z(x_0,eul_ang1)
            x_0 = rotate_X(x_0,eul_ang2)
            x_0 = rotate_Z(x_0,eul_ang3)
            y_0 = pylab.array([0.0,1.0,0.0])
            y_0 = rotate_Z(y_0,eul_ang1)
            y_0 = rotate_X(y_0,eul_ang2)
            y_0 = rotate_Z(y_0,eul_ang3)
            z_0 = pylab.array([1.0,0.0,0.0])
            z_0 = rotate_Z(z_0,eul_ang1)
            z_0 = rotate_X(z_0,eul_ang2)
            z_0 = rotate_Z(z_0,eul_ang3)
            N = self.map3d.shape[0]
            extra_space = 10
            map2d = pylab.zeros(shape=(N,N),dtype="complex64")
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
                        vector = pylab.array([(N-1)/2.0,(N-1)/2.0,(N-1)/2.0])-(N-1)/2.0*(x_0+y_0+z_0)+(z*z_0+y*y_0+x*x_0)
                        cases = [pylab.floor,lambda x: 1.0 + pylab.floor(x)]
                        for y_func in cases:
                            for x_func in cases:
                                for z_func in cases:
                                    if valid_coordinate(pylab.array([z_func(vector[0]),y_func(vector[1]),x_func(vector[2])])):
                                    #print "%i %i %i %i" % (z_func(vector[0]),y_func(vector[1]),x_func(vector[2]),N)
                                        map2d[y,x] +=\
                                            self.dX*\
                                            (1.0-abs(z_func(vector[0]) - vector[0]))*\
                                            (1.0-abs(y_func(vector[1]) - vector[1]))*\
                                            (1.0-abs(x_func(vector[2]) - vector[2]))*\
                                            self.map3d[int(z_func(vector[0])),
                                                  int(y_func(vector[1])),
                                                  int(x_func(vector[2]))]
        if "padding" in kwargs.keys():
            if kwargs["padding"]=="squared":
                N_new = max([map2d.shape[0],map2d.shape[1]])
                map2d_new = pylab.zeros(shape=(N_new,N_new))
                map2d_new[(N_new-map2d.shape[0])/2:(N_new-map2d.shape[0])/2+map2d.shape[0],
                          (N_new-map2d.shape[1])/2:(N_new-map2d.shape[1])/2+map2d.shape[1]] = map2d[:,:]
                map2d = map2d_new
        return map2d

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
            map3d = f.create_dataset('map3d', self.map3d.shape, self.map3d.dtype)
            map3d[:,:,:] = self.map3d[:,:,:]
            f['voxel_dimensions_in_m'] = self.dX
            f.close()
        else:
            print 'ERROR: Invalid filename extension, has to be \'.h5\'.'

    def load_map3d(self,filename):
        """
        Function loads refractive index map from an hdf5 file:
        ======================================================
        
        Arguments:
        
        - filename: Filename.

        """

        if filename[-3:] == '.h5':
            import h5py
            f = h5py.File(filename,'r')
            self.map3d = f['map3d'].value.copy()
            if f['voxel_dimensions_in_m'].value != self.dX: config.OUT.write("WARNING: Sampling of map and setup does not match.")
            self.dX = f['voxel_dimensions_in_m'].value
            f.close()
        else:
            print 'ERROR: Invalid filename extension, has to be \'.h5\'.'

    def get_area(self):
        """
        Function returns projected area of the sample along current beam axis:
        ======================================================================

        """
        projection = self.project(self.euler_angle_0,self.euler_angle_1,self.euler_angle_2)
        binary = pylab.ones_like(projection)
        binary[abs(projection)<pylab.median(abs(projection))] = 0
        area = binary.sum()*self.dX**2
        return area
        
class SampleSphere:
    """
    A class of the input-object.
    Sample is a homogeneous sphere defined by a radius and a material object.

    """

    def __init__(self,**kwargs):
        self.radius = kwargs.get('radius',225E-09)
        self._parent = kwargs.get('parent',None)

        material_kwargs = kwargs.copy()
        try: material_kwargs.pop('parent')
        except: pass
        try: material_kwargs.pop('radius')
        except: pass
        material_kwargs['parent'] = self
        self.material = Material(**material_kwargs)

    def get_area(self):
        """ Calculates area of projected sphere """
        return pylab.pi*self.radius**2
