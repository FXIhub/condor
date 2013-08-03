import pylab
import source, sample, detector, propagation, config, proptools



class Input:
    """
    Input object which is the necessary argument of function 'propagator'
    It contains all the information about the experimental setup (objects 'source', 'sample', 'detector')

    """
    
    def __init__(self,configfile=None):
        """
        Function initializes input-object:
        ==================================
        Arguments:
        - configfile: Filename of configuration file. If not given variables are set to default values.

        """
        self.source = source.Source(parent=self)
        self.sample = sample.SampleSphere(parent=self)
        self.detector = detector.Detector(parent=self)
        self.propagation = propagation.Propagation(parent=self)
        if configfile != None:
            self.read_configfile(configfile)
            config.OUT.write("Set configuration in accordance to given configuration-file: %s\n" % configfile)
        else:
            config.OUT.write("Initial values set to default values.\n")
    
    def set_sample_empty_map(self):
        """
        Function creates empty densitymap. Densitymap resolution is set according to the detector geometry.
        """
        self.sample = SampleMap(parent=self)

    def load_sample_map(self,filename,edge_length,**materialargs):
        """
        Creates empty densitymap. Size in accordance to given detector geometry.
        Densitymap resolution is set according to the detector geometry.
        """
        self.sample = SampleMap(parent=self)
        self.sample.load_map3d(filename)
        self.sample.dX = edge_length/(1.*(self.sample.map3d.shape[0]-1))
        if materialargs != {}:
            materialargs['parent'] = self.sample
            self.sample.material = Material(**materialargs)
            self.sample.map3d *= 1-self.sample.material.get_n()

    def set_sample_icosahedral_virus_map(self,radius=None,**kwargs):
        """
        Creates refractive index map of icosahedron. More information can be found in _makedm_icosahedron(...).
        """
        kwargs['materialtype'] = 'virus'
        self.set_sample_icosahedral_map(radius,**kwargs)

    def set_sample_icosahedral_map(self,radius=None,**kwargs):
        """
        Creates refractive index map of icosahedron. More information can be found in _makedm_icosahedron(...).
        """
        if radius == None:
            radius = self.sample.radius
        self.sample = SampleMap(parent=self)
        self.sample.put_icosahedral_sample(radius,
                                           **kwargs)
        self.sample.radius = radius

    def set_sample_spheroid_map(self,a,c,**kwargs):
        """
        Creates refractive index map of spheroid.
        """
        self.sample = SampleMap(parent=self)
        self.sample.put_spheroid_sample(a,
                                        c,
                                        **kwargs)
        self.sample.radius = (a**2*c)**(1/3.)

    def set_sample_sphere_map(self,radius=225E-09,**materialargs):
        """
        Creates refractive index map of sphere of given radius and material.
        """
        self.sample = SampleMap(parent=self)
        self.sample.put_sphere(radius,**materialargs)
        self.sample.radius = radius 

    def set_sample_homogeneous_sphere(self,**kwargs):
        """
        Sets sample object to homogeneous sphere having the given radius. Atomic composition values and massdensity are set according to the given material arguments.
        Examples for usage:
        - setting atomic composition values and massdensity manually:
          set_sample_homogeneous_sphere(diameter=100E-09,massdensity=1000,cH=2,cO=1)
        - setting atomic composition values and massdensity according to given materialtype:
          set_sample_homogeneous_sphere(diameter=100E-09,materialtype='protein')
        available materialtypes: 'protein', 'virus', 'cell', 'latexball', 'water'
        """ 
        self.sample = SampleSphere(parent=self,**kwargs)

    def set_sample_homogeneous_spheroid(self,**kwargs):
        self.sample = SampleSpheroid(parent=self,**kwargs)

    def read_configfile(self,configfile):
        """ 
        Function reads given configuration file and (over-)writes settings in the input-object.
        =======================================================================================
        
        Arguments:
        
        - configfile: Filename of the configuration file. [no default value]

        """
        C = ConfigParser.ConfigParser()
        try:
            C.readfp(open(configfile))
        except IOError:
            print "ERROR: Can't read configuration-file."
            return
        self.source.photon.set_wavelength(C.getfloat('source','wavelength'))
        self.source.focus_diameter = C.getfloat('source','focus_diameter')
        self.source.pulse_energy = C.getfloat('source','pulse_energy')

        args = {}
        args['distance'] = C.getfloat('detector','distance')
        args['pixel_size'] = C.getfloat('detector','pixel_size')
        args['binning'] = C.getint('detector','binning')
        args['Nx'] = C.getint('detector','Nx')
        args['Ny'] = C.getint('detector','Ny')
        if C.get('detector','cx') != 'middle': args['cx'] = C.getfloat('detector','cx')
        if C.get('detector','cy') != 'middle': args['cy'] = C.getfloat('detector','cy')
        args['saturation_level'] = C.getfloat('detector','saturation_level')
        if C.get('detector','mask') == 'none':
            args['x_gap_size_in_pixel'] = C.getint('detector','x_gap_size_in_pixel')
            args['y_gap_size_in_pixel'] = C.getint('detector','y_gap_size_in_pixel')
            args['hole_diameter_in_pixel'] = C.getint('detector','hole_diameter_in_pixel')
        else:
            M = pylab.imread(C.get('detector','mask'))
            M = M[:,:,0]
            M[abs(M)>0] = 1
            M[M==0] = 0
            args['mask'] = M.copy()
        self.detector = Detector(**args)
        self.propagation.rs_oversampling = C.getfloat('propagation','rs_oversampling')
        self.propagation.N_processes = C.getint('propagation','N_processes')

        sample_type = C.get('sample','sample_type','uniform_sphere')

        mat = C.get('sample','material_type','none')
        matargs = []
        if mat == 'none':
            pass
        elif mat == 'custom':
            if mat == 'custom':
                cX_list = C.items('sample')
                for cX_pair in cX_list:
                    if cX_pair[0][0] == 'c':
                        el = cX_pair[0]
                        el = el[1:].capitalize()
                        val = float(cX_pair[1])
                        matargs.append(("c%s" % el,val))
            matargs.append(('massdensity',C.getfloat('sample','massdensity')))
        else:
            matargs.append(('materialtype',mat))
        matargs= dict(matargs)
        
        
        if sample_type == 'uniform_sphere':
            matargs['diameter']=C.getfloat('sample','size')
            self.set_sample_homogeneous_sphere(**matargs)

        elif sample_type == 'uniform_spheroid':
            matargs['a']=C.getfloat('sample','a_diameter')/2.
            matargs['c']=C.getfloat('sample','c_diameter')/2.
            matargs['theta'] = C.getfloat('sample','theta')
            matargs['phi'] = C.getfloat('sample','phi')
            #euler_angle_2 = C.getfloat('sample','phi')
            self.set_sample_homogeneous_spheroid(**matargs)

        elif sample_type == 'map3d':
            geometry = C.get('sample','geometry','none')
            if geometry == 'none':
                size = C.getfloat('sample','size')
                self.set_sample_empty_map()
            elif geometry == 'icosahedron':
                size = C.getfloat('sample','size')
                euler_angle_0 = C.getfloat('sample','theta')
                euler_angle_1 = C.getfloat('sample','psi')
                euler_angle_2 = C.getfloat('sample','phi')
                self.set_sample_icosahedral_map(radius=size/2.,euler_angle_0=euler_angle_0,euler_angle_1=euler_angle_1,euler_angle_2=euler_angle_2,**matargs)
            elif geometry == 'spheroid':
                a = C.getfloat('sample','a_diameter')/2.
                c = C.getfloat('sample','c_diameter')/2.
                euler_angle_0 = C.getfloat('sample','theta')
                euler_angle_1 = C.getfloat('sample','psi')
                euler_angle_2 = C.getfloat('sample','phi')
                self.set_sample_spheroid_map(a,c,euler_angle_0=euler_angle_0,euler_angle_1=euler_angle_1,euler_angle_2=euler_angle_2,**matargs)
            else:
                size = C.getfloat('sample','size')
                if mat == 'none':
                    self.load_sample_map(geometry,size)
                else:
                    self.load_sample_map(geometry,size,**matargs)
            
            

    def get_real_space_resolution_element(self):
        """
        Function returns real space sampling.
        =====================================
        Formula: dX = pi / qmax

        """
        dX = pylab.pi / self.get_absq_max()
        return dX

    def get_max_achievable_crystallographic_resolution(self):
        """
        Function returns maximal horizontal and vertical crystallographoc resolution achievable in the current setup.
        =============================================================================================================
        Formula: Rx_c = dx * 2 = 2 * pi / qxmax ; Ry_c = dy * 2 = 2 * pi / qymax

        """
        dx = 2 * pylab.pi / self.get_absqx_max()
        dy = 2 * pylab.pi / self.get_absqy_max()
        return [dx,dy]

    def generate_absqmap(self):
        X,Y = pylab.meshgrid(pylab.arange(self.detector.mask.shape[1]),
                             pylab.arange(self.detector.mask.shape[0]))
        # THIS CAST IS VERY IMPORTANT, in python A += B is not the same as A = A + B
        X = pylab.float64(X)
        Y = pylab.float64(Y)
        X -= self.detector.get_cx('binned')
        Y -= self.detector.get_cy('binned')
        p = self.detector.get_pixel_size('binned')
        D = self.detector.distance
        w = self.source.photon.get_wavelength()
        return proptools.generate_absqmap(X,Y,p,D,w)

    def generate_qmap(self,nfft_scaled=False):
        X,Y = pylab.meshgrid(pylab.arange(self.detector.mask.shape[1]),
                             pylab.arange(self.detector.mask.shape[0]))
        # THIS CAST IS VERY IMPORTANT, in python A += B is not the same as A = A + B
        X = pylab.float64(X)
        Y = pylab.float64(Y)
        X -= self.detector.get_cx('binned')
        Y -= self.detector.get_cy('binned')
        p = self.detector.get_pixel_size('binned')
        D = self.detector.distance
        w = self.source.photon.get_wavelength()
        try: E0 = self.euler_angle_0
        except: E0 = 0.
        try: E1 = self.euler_angle_1
        except: E1 = 0.
        try: E2 = self.euler_angle_2
        except: E2 = 0.
        qmap = proptools.generate_qmap(X,Y,p,D,w,E0,E1,E2)
        if nfft_scaled == True:
            return qmap/self.get_absq_max()*0.5
        else:
            return qmap
        
    def get_phase_ramp(self,qmap,dvector):
        return pylab.exp(-1.j*(dvector[0]*qmap[:,:,0]+
                               dvector[1]*qmap[:,:,1]+
                               dvector[2]*qmap[:,:,2]))

    def get_absqx_max(self):
        wavelength = self.source.photon.get_wavelength()
        x_max = max([self.detector.get_cx('binned'),self.detector.mask.shape[1]-1-self.detector.get_cx('binned')]) * self.detector.get_pixel_size('binned')
        R_Ewald = 2*pylab.pi/wavelength
        phi = pylab.arctan2(x_max,self.detector.distance)
        return 2 * R_Ewald * pylab.sin(phi/2.0)

    def get_absqy_max(self):
        wavelength = self.source.photon.get_wavelength()
        y_max = max([self.detector.get_cy('binned'),self.detector.mask.shape[1]-1-self.detector.get_cy('binned')]) * self.detector.get_pixel_size('binned')
        R_Ewald = 2*pylab.pi/wavelength
        phi = pylab.arctan2(y_max,self.detector.distance)
        return 2 * R_Ewald * pylab.sin(phi/2.0)

    def get_absqz_max(self):
        absqx_max = self.get_absqx_max()
        absqy_max = self.get_absqy_max()
        wavelength = self.source.photon.get_wavelength()
        R_Ewald = 2*pylab.pi/wavelength
        phi = pylab.arcsin(pylab.sqrt(absqx_max**2+absqy_max**2)/R_Ewald)
        return R_Ewald * (1-pylab.cos(phi))

    def get_absq_max(self):
        #return pylab.sqrt(self.get_absqx_max()**2+self.get_absqy_max()**2+self.get_absqz_max()**2)
        return max([self.get_absqx_max(),self.get_absqy_max()])
