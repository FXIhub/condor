import sys,pylab,time,multiprocessing
if "utils" not in sys.path: sys.path.append("utils")
import config,imgutils,proptools

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
            self.map3d = pylab.zeros(shape=(N,N,N),dtype="complex64")
        self._parent = kwargs.get('parent',None)
        self.radius = kwargs.get('radius',None)
        self.dX = kwargs.get('dX',self._parent.get_real_space_resolution_element()/(1.0*self._parent.propagation.rs_oversampling))
        self.euler_angle_0 = kwargs.get('euler_angle_0',0.0)
        self.euler_angle_1 = kwargs.get('euler_angle_1',0.0)
        self.euler_angle_2 = kwargs.get('euler_angle_2',0.0)

    def set_random_orientation(self):
        [e0,e1,e2] = proptools.random_euler_angles()
        self.euler_angle_0 = e0
        self.euler_angle_1 = e1
        self.euler_angle_2 = e2
        
    def put_custom_map(self,map_add,x=0.,y=0.,z=0.,mode='refill'):
        Nx = self.map3d.shape[2]
        Ny = self.map3d.shape[1]
        Nz = self.map3d.shape[0]
        origin = pylab.array([(self.map3d.shape[0]-1)/2.0,
                              (self.map3d.shape[1]-1)/2.0,
                              (self.map3d.shape[2]-1)/2.0])
        zmin = round(origin[0]+z-map_add.shape[0]/2.0)
        ymin = round(origin[1]+y-map_add.shape[1]/2.0)
        xmin = round(origin[2]+x-map_add.shape[2]/2.0)
        zoff = 0
        yoff = 0
        xoff = 0
        if zmin < 0:
            Nz += -zmin
            zoff = -zmin
        if ymin < 0:
            Ny += -ymin
            yoff = -ymin
        if xmin < 0:
            Nx += -xmin
            xoff = -xmin
        change_map_size = False
        if self.map3d.shape[0] < zmin+map_add.shape[0]:
            Nz += (zmin+map_add.shape[0])-self.map3d.shape[0]
            change_map_size = True
            print Nz
        if self.map3d.shape[1] < ymin+map_add.shape[1]:
            Ny += (ymin+map_add.shape[1])-self.map3d.shape[1]
            change_map_size = True
            print Ny
        if self.map3d.shape[2] < xmin+map_add.shape[2]:
            Nx += (xmin+map_add.shape[2])-self.map3d.shape[2]
            change_map_size = True
            print Nx
        if change_map_size:
            temp = self.map3d.copy()
            self.map3d = pylab.zeros(shape=(Nz,Ny,Nx),dtype="complex64")
            self.map3d[zoff:zoff+temp.shape[0],yoff:yoff+temp.shape[1],xoff:xoff+temp.shape[2]] = temp[:,:,:]
        #self.map3d[zoff+zmin:zoff+zmin+map_add.shape[0],
        #           yoff+ymin:yoff+ymin+map_add.shape[1],
        #           xoff+xmin:xoff+xmin+map_add.shape[2]] += map_add[:,:,:]
        if mode=='fill':
            M = self.map3d[zoff+zmin:zoff+zmin+map_add.shape[0],
                           yoff+ymin:yoff+ymin+map_add.shape[1],
                           xoff+xmin:xoff+xmin+map_add.shape[2]]
            M[M==0.] = map_add[M==0.]
            M[(M<M.max())*(M!=0.)] = (map_add[(M<M.max())*(M!=0.)]+M[(M<M.max())*(M!=0.)])/2.
        elif mode=='refill':
            self.map3d[zoff+zmin:zoff+zmin+map_add.shape[0],
                       yoff+ymin:yoff+ymin+map_add.shape[1],
                       xoff+xmin:xoff+xmin+map_add.shape[2]] = map_add[:,:,:]
        elif mode=='add':
            self.map3d[zoff+zmin:zoff+zmin+map_add.shape[0],
                       yoff+ymin:yoff+ymin+map_add.shape[1],
                       xoff+xmin:xoff+xmin+map_add.shape[2]] += map_add[:,:,:]


    # TESTING REQUIRED
    def smooth_map3d(self,factor):
        self.map3d = imgutils.smooth3d(self.map3d,factor)

    def put_spheres_stack(self,Nppp,width,height,sphere_diameter,layer_thickness,**materialargs):
        material_obj = Material(self,**materialargs)
        material_obj_water = Material(self,materialtype='water')
        dn = 1.0-material_obj.get_n()
        dn_water = 1.0-material_obj_water.get_n()
        edge_pix = width/self.dX
        diameter_pix  = sphere_diameter/self.dX
        thickness_pix = layer_thickness/self.dX
        Np = int(pylab.ceil(height/layer_thickness))
        map3d = imgutils.generate_random_colloid_planes(edge_pix,diameter_pix,thickness_pix,Np,Nppp) * dn
        map3d[abs(map3d)<abs(dn_water)] = dn_water
        if self.map3d.max() == 0:
            self.map3d = map3d
        else:
            self.put_custom_map(map3d)
    
    def put_into_water(self):
        M = Material(self,materialtype='water')
        dn = 1.-M.get_n()
        self.map3d.real += dn.real*(1-self.map3d.real/self.map3d.real.max())
        self.map3d.imag += dn.imag*(1-self.map3d.imag/self.map3d.imag.max())
        self.map3d -= dn

    def put_periodic_filament(self,diameter,layer_thickness,length,**materialargs):
        Nx = length/self.dX
        Ny = diameter/self.dX
        Nz = diameter/self.dX
        T = layer_thickness/self.dX
        Z,Y,X = 1.0*pylab.mgrid[0:Nz,0:Ny,0:Nx]
        Z -= Nz/2.
        X -= Nx/2.
        Y -= Ny/2.
        R = pylab.sqrt(Y**2+Z**2)
        map3d = R.copy()
        map3d[R<Ny/2.] = 1-R[R<Ny/2.]/(Ny/2.)
        map3d[R>=Ny/2.] = 0
        map3d *= pylab.sin(2*pylab.pi*X/T)
        material_obj = Material(self,**materialargs)
        dn = 1.0-material_obj.get_n()
        self.map3d = dn*map3d

    def put_sphere(self,radius,x=0,y=0,z=0,fillmode='refill',**materialargs):
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
        spheremap = pylab.zeros(shape=R.shape,dtype="complex64")
        spheremap[R<R_N] = 1
        spheremap[abs(R_N-R)<0.5] = 0.5+0.5*(R_N-R[abs(R_N-R)<0.5])
        spheremap *= dn
        if self.map3d.max() == 0 and not x and not y and not z:
            self.map3d = spheremap
        else:
            x_N = x/self.dX
            y_N = y/self.dX
            z_N = z/self.dX
            self.put_custom_map(spheremap,x_N,y_N,z_N,fillmode)
 
    def put_spheroid(self,a,b,x=None,y=None,z=None,eul_ang0=0.0,eul_ang1=0.0,eul_ang2=0.0,**materialargs):
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
        config.OUT.write("Rotate grid.\n")
        [X,Y,Z] = imgutils.rotate_3d_grid(X,Y,Z,eul_ang0,eul_ang1,eul_ang2)
        config.OUT.write("Build sphorid on a grid of %i x %i x %i\n" % (size,size,size))
        spheroidmap = (X**2+Y**2)/a_N**2+Z**2/b_N**2
        spheroidmap[spheroidmap<=1] = 1
        spheroidmap[spheroidmap>1] = 0
        smooth_factor = 1.5
        config.OUT.write("Smoothing by a factor of %f\n" % smooth_factor)
        spheroidmap = abs(imgutils.smooth3d(spheroidmap,smooth_factor))
        spheroidmap /= spheroidmap.max()
        spheroidmap = pylab.array(spheroidmap,dtype="complex64")
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
        self.put_sphere(radius,x,y,z,'fill',cAu=1,massdensity=config.DICT_massdensity['Au'])        

    def put_icosahedral_virus(self,radius,x=0.,y=0.,z=0.,**kwargs):
        kwargs_cp = kwargs.copy()
        kwargs_cp['materialtype'] = 'virus'
        dn = self._makedm_icosahedron(radius,**kwargs_cp)
        if self.map3d.max() == 0.0 and x==0. and y==0. and z==0.:
            self.map3d = dn
        else:
            x_N = x/self.dX
            y_N = y/self.dX
            z_N = z/self.dX
            self.put_custom_map(dn,x_N,y_N,z_N)

    def put_icosahedral_sample(self,radius,x=0.,y=0.,z=0.,**kwargs):
        kwargs_cp = kwargs.copy()
        kwargs_cp['materialtype'] = kwargs['materialtype']
        dn = self._makedm_icosahedron(radius,**kwargs_cp)
        if self.map3d.max() == 0.0 and x==0. and y==0. and z==0.:
            self.map3d = dn
        else:
            x_N = x/self.dX
            y_N = y/self.dX
            z_N = z/self.dX
            self.put_custom_map(dn,x_N,y_N,z_N)

    def make_grid(self,Nx,Ny,Nz,spacing=2,thickness=1):
        G = pylab.ones(shape=(Nz,Ny,Nx))
        for iz in pylab.arange(0,Nz,spacing):
            for it in pylab.arange(0,thickness,1):
                G[iz+it,:,:] = 0
        return G

    def plot_map3d(self):
        from enthought.mayavi import mlab
        s = mlab.pipeline.scalar_field(abs(self.map3d))
        plane_1 = mlab.pipeline.image_plane_widget(s,plane_orientation='x_axes',
                                                   slice_index=self.map3d.shape[2]/2)
        plane_2 = mlab.pipeline.image_plane_widget(s,plane_orientation='y_axes',
                                                   slice_index=self.map3d.shape[1]/2)
        mlab.show()

    def plot_fmap3d(self):
        from enthought.mayavi import mlab
        import spimage
        M = spimage.sp_image_alloc(self.map3d.shape[0],self.map3d.shape[1],self.map3d.shape[2])
        M.image[:,:,:] = self.map3d[:,:,:]
        M.mask[:,:,:] = 1
        fM = spimage.sp_image_fftw3(M)
        fM.mask[:,:,:] = 1
        fsM = spimage.sp_image_shift(fM)
        self.fmap3d = abs(fsM.image).copy()
        self.fmap3d[self.fmap3d!=0] = pylab.log10(self.fmap3d[self.fmap3d!=0])
        
        print pylab.log10(abs(self.fmap3d))
        s = mlab.pipeline.scalar_field(self.fmap3d)
        plane_1 = mlab.pipeline.image_plane_widget(s,plane_orientation='x_axes',
                                                   slice_index=self.fmap3d.shape[2]/2)
        plane_2 = mlab.pipeline.image_plane_widget(s,plane_orientation='y_axes',
                                                   slice_index=self.fmap3d.shape[1]/2)
        mlab.show()


    def _makedm_icosahedron(self,radius,**kwargs):  
        """
        Function returns a refractive index map of a homogeneous icosahedron:
        =====================================================================

        Arguments:

        - radius: Radius in meter (volume equals volume of a sphere with given radius). [no default value]
        
        Keyword arguments:

        - euler1,euler2,euler3: orientation defined by Euler angles in rad. [0.0,0.0,0.0]
        - Nlayers: generate icosahedron that is radially layered with Nlayers across.
          If 'Nlayers' is set to 0 icosahedron is being homogeneously filled. [0]

        EITHER:
        - materialtype: Predefined materialtype ('protein','virus','cell','latexball','water',... For a complete list look up \'config.py\'.)
        
        OR:
        - c\'X\': fraction of element \'X\' of atomic number density in material (normalization is done internally).
        - massdensity: Massdensity of material in kg/m^3

        """

        if 'massdensity' in kwargs:
            materialargs = {'massdensity':kwargs['massdensity']}
            for key in kwargs.keys():
                if kwargs[i][0] == 'c': materialargs[key] = kwargs[key]
        else:
            materialargs = {'materialtype':kwargs['materialtype']}
        material_obj = Material(self,**materialargs)
        dn = 1.0 - material_obj.get_n()

        euler1 = kwargs.get('euler1',0.0)
        euler2 = kwargs.get('euler2',0.0)
        euler3 = kwargs.get('euler3',0.0)
        Nlayers = kwargs.get('Nlayers',0)

        t_start = time.time()

        a = radius*(16*pylab.pi/5.0/(3+pylab.sqrt(5)))**(1/3.0)
        Rmax = pylab.sqrt(10.0+2*pylab.sqrt(5))*a/4.0 # radius at corners
        Rmin = pylab.sqrt(3)/12*(3.0+pylab.sqrt(5))*a # radius at faces
        nRmax = Rmax/self.dX
        nRmin = Rmin/self.dX

        # smooth
        s = 1.0

        N = int(pylab.ceil(2.2*(nRmax)))
        
        config.OUT.write("... build icosahedral geometry ...\n")
        
        n_list = imgutils.get_icosahedron_normal_vectors(euler1,euler2,euler3)
        X,Y,Z = 1.0*pylab.mgrid[0:N,0:N,0:N]
        X -= (N-1)/2.
        Y -= (N-1)/2.
        Z -= (N-1)/2.

        config.OUT.write("... %i x %i x %i grid (%i voxels) ...\n" % (N,N,N,N**3))
        icomap = pylab.zeros((len(n_list),N,N,N))
        
        if Nlayers == 0:

            for i in range(len(n_list)):
                icomap[i,:,:,:] = -(X*n_list[i][2]+Y*n_list[i][1]+Z*n_list[i][0]-nRmin)

            M = icomap.copy()
            #icomap[abs(M)<=0.5*s] = abs(0.5+icomap[abs(M)<=0.5*s]/s)
            icomap[M<=0.5*s] = 0
            icomap[M>0.5*s] = 1
            icomap = icomap.min(0)

        else:
            nlayerdistance = 2*nRmin/(1.0*Nlayers)

            for i in range(len(n_list)):
                icomap[i,:,:,:] = X*n_list[i][2]+Y*n_list[i][1]+Z*n_list[i][0]

            M = icomap-nRmin
            M = M.max(0)
            M[M>0] = 0

            icomap = pylab.zeros_like(M)
            icomap = (1-pylab.cos(2*pylab.pi*M/nlayerdistance))/2.
            #f = 1.
            #icomap[abs(M)<nRmin*f] = (1-pylab.cos(2*pylab.pi*M[abs(M)<nRmin*f]/nlayerdistance))/2.
            #icomap[abs(M)>=nRmin*3/4.] = (1-pylab.cos(2*pylab.pi*Z[abs(M)>=nRmin*3/4.]/nlayerdistance))/2.

        t_stop = time.time()

        config.OUT.write("Time: %f sec\n" % (t_stop-t_start))

        return dn*icomap    

    def project(self,eul_ang0=None,eul_ang1=None,eul_ang2=None,**kwargs):
        """ Projection of 3-dimensional map."""
        if not eul_ang0: eul_ang0 = self.euler_angle_0
        if not eul_ang1: eul_ang1 = self.euler_angle_1
        if not eul_ang2: eul_ang2 = self.euler_angle_2

        if eul_ang0 == 0.0 and eul_ang1 == 0.0 and eul_ang2 == 0.0:
            map2d = pylab.zeros((self.map3d.shape[1],self.map3d.shape[2]))
            for iy in pylab.arange(0,map2d.shape[0]):
                for ix in pylab.arange(0,map2d.shape[1]):
                    map2d[iy,ix] = self.map3d[:,iy,ix].real.sum()*self.dX*pylab.exp(-self.map3d[:,iy,ix].imag.sum()*self.dX)
        else:
            x_0 = proptools.rotation(pylab.array([0.0,0.0,1.0]),eul_ang0,eul_ang1,eul_ang2)
            y_0 = proptools.rotation(pylab.array([0.0,1.0,0.0]),eul_ang0,eul_ang1,eul_ang2)
            z_0 = proptools.rotation(pylab.array([1.0,0.0,0.0]),eul_ang0,eul_ang1,eul_ang2)
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

    def get_area(self,mode='auto'):
        """
        Function returns projected area of the sample along current beam axis:
        ======================================================================

        """
        if mode == 'auto':
            if self.radius != None: return pylab.pi*self.radius**2
        projection = self.project(self.euler_angle_0,self.euler_angle_1,self.euler_angle_2)
        binary = pylab.ones(shape=projection.shape)
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


class SampleSpheres:
    """
    A class of the input-object.
    Sample is a set of homogeneous spheres defined by a radii and material objects.

    """

    def __init__(self,**kwargs):
        self._parent = kwargs.get('parent',None)

        material_kwargs = kwargs.copy()
        try: material_kwargs.pop('parent')
        except: pass
        material_kwargs['parent'] = self
        if len(material_kwargs) != 0: self.material = Material(**material_kwargs)
        else: self.material = None

        self.x = []
        self.y = []
        self.z = []
        self.r = []
        self.m = []

    def fill_cube_randomly(self,a,N,d):
        X = imgutils.get_random_circle_positions(N,d/a,3)
        z = X[:,0]
        y = X[:,1]
        x = X[:,2]
        for i in range(N):
            self.add_sphere(d/2.,a*x[i],a*y[i],a*z[i])

    def fill_icosahedron_randomly(self,r,N,d):
        Ni = 0
        ns = imgutils.get_icosahedron_normal_vectors()
        f = 1.2
        while Ni < N:
            to_add_indices = []
            Nnow = int(N*f)
            X = imgutils.get_random_circle_positions(Nnow,d/r/2.,3)
            X = (X-0.5)*2.
            z = X[:,0]
            y = X[:,1]
            x = X[:,2]
            for i in range(Nnow):
                outside = False
                for n in ns:
                    if pylab.dot(X[i,:],n) > 0.7:
                        outside = True
                        break
                if not outside:
                    to_add_indices.append(i)
                    Ni += 1
                    print Ni
            f += 0.2
        for i in to_add_indices:
            self.add_sphere(d/2.,r*x[i],r*y[i],r*z[i])

    def fill_icosahedron_radial_chains(self,r,N,d,t,s):
        ns = imgutils.get_icosahedron_normal_vectors()
      
        def pick_random_point():
            [u,v] = liSst(pylab.rand(2))
            theta = 2*pylab.pi*u
            phi = pylab.arccos(2*v-1)
            
            
            return distances.min()

        X = pylab.zeros(shape=(N,3))

        D = d/2./r
        Dsq = D**2

        i = 0
        intersect = False
        while i < N:
            X[i,:] = pylab.rand(3)[:]
            X[i] = (X[i]-0.5)*2.
            for j in range(i):
                intersect = False
                if pylab.dot(X[i]-X[j],X[i]-X[j]) <= Dsq:
                    intersect = True
                    print intersect
                    break
            if intersect == False:
                for n in ns:
                    outside = False
                    if pylab.dot(X[i,:],n) > 0.7:
                        outside = True
                        break
                if not outside:
                    if abs((get_distance(r*X[i]) % t)-t/2.) < s:
                        self.add_sphere(d/2.,r*X[i,2],r*X[i,1],r*X[i,0])
                        i += 1
                        print i

    def fill_icosahedron_radial_layers(self,r,N,d,t,s):
        ns = imgutils.get_icosahedron_normal_vectors()
      
        def get_distance(pos):
            distances = pylab.zeros(len(ns))
            for i in range(len(ns)):
                distances[i] = pylab.dot(ns[i],pos)-pylab.sqrt(pylab.dot(ns[i],ns[i]))
            return distances.min()

        X = pylab.zeros(shape=(N,3))

        D = d/2./r
        Dsq = D**2

        i = 0
        intersect = False
        while i < N:
            X[i,:] = pylab.rand(3)[:]
            X[i] = (X[i]-0.5)*2.
            for j in range(i):
                intersect = False
                if pylab.dot(X[i]-X[j],X[i]-X[j]) <= Dsq:
                    intersect = True
                    print intersect
                    break
            if intersect == False:
                for n in ns:
                    outside = False
                    if pylab.dot(X[i,:],n) > 0.7:
                        outside = True
                        break
                if not outside:
                    if abs((get_distance(r*X[i]) % t)-t/2.) < s:
                        self.add_sphere(d/2.,r*X[i,2],r*X[i,1],r*X[i,0])
                        i += 1
                        print i
              
    def fill_icosahedron_parallel_chains(self,r,N,d,t,axis=2):
        ns = imgutils.get_icosahedron_normal_vectors()
      
        def get_distance(pos):
            distances = pylab.zeros(len(ns))
            for i in range(len(ns)):
                distances[i] = pylab.dot(ns[i],pos)-pylab.sqrt(pylab.dot(ns[i],ns[i]))
            return distances.min()

        X = pylab.zeros(shape=(N,3))

        D = d/(2.*r)
        Dsq = D**2

        i = 0
        intersect = False
        while i < N:
            X[i,:] = pylab.rand(3)[:]
            X[i] = (X[i]-0.5)*2.

            for j in range(i):
                intersect = False
                if pylab.dot(X[i]-X[j],X[i]-X[j]) <= Dsq:
                    intersect = True
                    print 'Missed insert'
                    break

            Xtemp = (X[i]-0.5)*2.
            istart = i

            outside_top = False
            while outside_top == False and i < N:
                if abs(Xtemp[0])>1. and abs(Xtemp[1])>1. and abs(Xtemp[2])>1.:
                    break
                else:
                    for n in ns:
                        if pylab.dot(Xtemp[:],n) > 0.7:
                            outside_top = True
                            break
                    if not outside_top:
                        self.add_sphere(d/2.,r*Xtemp[2],r*Xtemp[1],r*Xtemp[0])
                        i += 1
                        Xtemp[axis] += t/(2.*r)
                        print i

            if i < N:
                Xtemp = (X[istart]-0.5)*2.
                Xtemp[0] -= t/(2.*r)
            else: break

            outside_bot = False
            while outside_bot == False and i < N:
                if abs(Xtemp[0])>1. and abs(Xtemp[1])>1. and abs(Xtemp[2])>1.:
                    break
                for n in ns:
                    if pylab.dot(Xtemp[:],n) > 0.7:
                        outside_bot = True
                        break
                if not outside_bot:
                    self.add_sphere(d/2.,r*Xtemp[2],r*Xtemp[1],r*Xtemp[0])
                    i += 1
                    Xtemp[axis] -= t/(2.*r)
                    print i

    def make_filaments(self,diameter_element,N,length_filament,diameter_filament):
        X = pylab.zeros(shape=(N,2))

        i = 0
        intersect = False
        dsq = (diameter_element/diameter_filament)**2

        while i < N:
            X[i,:] = pylab.rand(2)[:]
            X[i,:] = (X[i,:]-0.5)*2.
            intersect = False
            for j in range(i):
                if pylab.dot(X[i,:]-X[j,:],X[i,:]-X[j,:]) <= dsq:
                    print 'Missed insert'
                    intersect = True
                    break
            if intersect == False:
                print int(round(length_filament/diameter_element))
                for j in range(int(round(length_filament/diameter_element))):
                    print 'adding sphere'
                    self.add_sphere(diameter_element/2.,
                                    diameter_element*j,
                                    X[i,0]*diameter_filament,
                                    X[i,1]*diameter_filament)
                i += 1
  
    def clear(self):
        self.x = []
        self.y = []
        self.z = []
        self.r = []
        self.m = []        

    def add_sphere(self,radius,x,y,z,material=None):
        self.x.append(x)
        self.y.append(y)
        self.z.append(z)
        self.r.append(radius)
        if material == None: material = self.material
        self.m.append(material)

    def plot_spheres(self):
        import mayavi.mlab
        mayavi.mlab.points3d(self.x,self.y,self.z,self.r,colormap="copper", scale_factor=1.)

    def get_area(self):
        return 1

    def get_radii(self):
        radii = []
        temp = pylab.zeros(len(self.r))
        while sum(temp) != len(temp):
            radii.append(pylab.array(self.r)[temp==0].min())
            temp[pylab.array(self.r)==radii[-1]] = 1
        return radii

    
