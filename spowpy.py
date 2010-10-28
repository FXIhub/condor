#! /usr/bin/python

# SPOWPY: Scattering experiment simulator for spheres and customized densitymaps
# Please type 'help spow()' for further information.

# Import libs
#------------
import pylab, sys, ConfigParser, numpy, types, pickle

ELEMENTS_FILE = open('elements.dat','r')
DICT_atomic_mass,DICT_scattering_factors = pickle.load(ELEMENTS_FILE)

F_MIN_ENERGY_EV = 0
F_MAX_ENERGY_EV = 0

for var in DICT_scattering_factors.values():
    if F_MIN_ENERGY_EV < var[0,0] or F_MIN_ENERGY_EV == 0: F_MIN_ENERGY_EV = var[0,0]
    if F_MAX_ENERGY_EV > var[-1,0] or F_MAX_ENERGY_EV == 0: F_MAX_ENERGY_EV = var[-1,0]

# Define global variables
#------------------------
MODE_SAMPLE_DENSITYMAP = 1
MODE_SAMPLE_HOMOGENEOUSSPHERE = 2
MODE_PRINT_NONE = 0
MODE_PRINT_STD = 1
MODE_PRINT_OBJ = 2
MODE_PLOT_NONE = 0
MODE_PLOT_1D = 1
MODE_PLOT_2D = 2
MODE_SCALING_PIX = 0
MODE_SCALING_NYQUISTPIX = 1

# Define global dictionaries
#---------------------------
# Typical realative atomic compositions (order: H,C,N.O,P,S,Au)
DICT_atomic_composition = {'protein':[86,52,13,15,0,3],'virus':[72.43,47.52,13.55,17.17,1.11,0.7],'cell':[23,3,1,10,0,1],'latexball':[1,1,0,0,0,0],'water':[2,0,0,1,0,0]}
# Typical realative atomic compositions (order: H,C,N.O,P,S)
DICT_massdensity = {'protein':1350,'virus':1455,'cell':1000,'latexball':1050,'water':998,'Au':19300}
# Physical constants [SI-units]
DICT_physical_constants = {'e':1.602176487E-19,'c':299792458,'h':6.62606896E-34,'re':2.8179402894E-15,'barn':1E-28,'u':1.66053886E-27}

# Input object in configuration
#------------------------------
class Input:
    """ INPUT CONFIGURATION OBJECT FOR 'spow()'\n\n
    An 'Input'-object is the (only) necessary argument for 'spow()'.\n\n
    The 'Input'-object contains all the information about:
    - experimental setup (in objects 'self.source', 'self.sample' and 'self.detector')
    - output- and program-settings of 'spow()' (in the object "it-'self'")\n\n
    Initialization of input-object:\n
    - <input_obj> = Input() 
      -> creates <input_obj> and sets all values to default values\n
    - <input_obj> = Inout(<conffilename>)
      -> creates <input_obj>. All variables are set to the values specified in the given configuration file"""
    def __init__(self,configfile=None,sample='homogeneous sphere',printmode='loglist',plotmode='none'):
        self.set_printmode(printmode)
        self.output = Output()
        self.output._tmploglist = _WritableObject()
        if printmode == 'stdout':
            clout = sys.stdout
        elif printmode == 'loglist':
            clout = self.output._tmploglist
        else:
            print "ERROR: %s is no valid printmode." % printmode
            return

        self.set_plotmode(plotmode)

        self.source = Source(self)

        if sample == 'homogeneous sphere':
            self.set_sample_homogeneous_sphere(self)
        elif sample == 'none':
            pass
        else:
            print "ERROR: %s is not a valid samplemode." % sample
            return

        self.detector = Detector(self)

        if configfile != None:
            self.read_configfile(configfile)
            clout.write("... set configuration in accordance to given configuration-file: %s ...\n" % configfile)
        else:
            clout.write("... initial values set to default values ...\n")
    
    def set_printmode(self,mode):
        """ Printmode determines the target of commandline outputs:
        mode: - 'stdout'  Comments will be printed to standard output
              - 'loglist' Comments will be printed to <output-object>.loglist"""
        if mode == 'stdout':
            self._printmode = MODE_PRINT_STD
        elif mode == 'loglist':
            self._printmode = MODE_PRINT_OBJ
        else:
            print "ERROR: No valid given argument."

    def set_plotmode(self,mode):
        """ Plotmode determines if graphs/patterns will be generated automatically after running 'spow()':
        mode: - 'none' Switches off plotting (default)
              - '1d'   1-dimensional graphs will be generated. They will show the distribution of scattered photons per pixel / Nyquist-pixel.
              - '2d'   2-dimensional patterns will be generated. They will show the distribution of scattered photons per pixel / Nyquist-pixel."""
        if mode == 'none':
            self._plotmode = MODE_PLOT_NONE
        elif mode == '1d':
            self._plotmode = MODE_PLOT_1D
        elif mode == '2d':
            self._plotmode = MODE_PLOT_2D
        else:
            print "ERROR: No valid argument."


    def set_sample_empty_densitymap(self,size):
        """
        Create empty densitymap. Edgelengths are specified by given argument 'size' in m.
        Densitymap resolution is set according to the detector geometry.
        """
        self.sample = DensitymapSample(self,size)

    def set_sample_virus_densitymap(self,radius,eul_ang1,eul_ang2,eul_ang3,speedup_factor=1):
        """
        Create virus of sphere-volume-equivalent given radius, rotates according to given Euler-angles euler_angle1, euler_angle2 and euler_angle3 [rad].
        Usage: set_sample_virus_densitymap(radius,euler-angle1,euler-angle2,euler-angle3,[speedup_factor])
        Densitymap resolution is set to highest resolution that can be achieved by the given detector geometry.
        For rough simulations that can be changed by setting the optional argument 'speedup_factor' to an integer bigger than 1.
        """
        self.sample = DensitymapSample(self,None)
        material_obj = Material(self.sample,materialtype='virus')
        [self.sample.densitymap2d,self.sample.densitymap3d,self.sample.densitymap_d] = self.sample._makedm_icosahedron(radius,material_obj,eul_ang1,eul_ang2,eul_ang3,speedup_factor)

    def set_sample_homogeneous_sphere(self,radius=225E-09,**materialargs):
        """
        Sets sample object to homogeneous sphere having the given radius. Atomic composition values and massdensity are set according to the given material arguments.
        Examples for usage:
        - setting atomic composition values and massdensity manually:
          set_sample_homogeneous_sphere(radius,massdensity=1000,cH=2,cO=1)
        - setting atomic composition values and massdensity according to given materialtype:
          set_sample_homogeneous_sphere(radius,materialtype='protein')
        available materialtypes: 'protein', 'virus', 'cell', 'latexball', 'water', 'Au'
        """ 
        self.sample = HomogeneoussphereSample(radius,**materialargs)

    def read_configfile(self,configfile):
        """ Reads given configuration file and sets configuration to the input-object """
        e = DICT_physical_constants['e']
        c = DICT_physical_constants['c']
        h = DICT_physical_constants['h']
        re = DICT_physical_constants['re']
        u = DICT_physical_constants['u']
        barn = DICT_physical_constants['barn']
        config = ConfigParser.ConfigParser()
        try:
            config.readfp(open(configfile))
        except IOError:
            print "ERROR: Can't read configuration-file."
            return
        self.source.wavelength = config.getfloat('source','wavelength')
        self.source.sizex = config.getfloat('source','sizex')
        self.source.sizey = config.getfloat('source','sizey')
        self.source.energy = config.getfloat('source','energy')
        mat = config.get('sample','material')
        args = []
        if mat == 'custom':
            cX_list = config.items('sample')
            for cX_pair in cX_list:
                if cX_pair[0][0] == 'c':
                    el = cX_pair[0]
                    el = el[1:].capitalize()
                    val = float(cX_pair[1])
                    args.append(("'c%s'" % el,val))
            args.append(('massdensity',config.getfloat('sample','massdensity')))
        else:
            keys = ['cH','cN','cO','cP','cS']
            for i in range(0,len(keys)):
                args.append((keys[i],DICT_atomic_composition[mat][i]))
            args.append(('massdensity',DICT_massdensity[mat]))
        args= dict(args)
        self.sample = HomogeneoussphereSample(self,config.getfloat('sample','radius'),**args)
        self.detector.distance = config.getfloat('detector','distance')
        self.detector.psize = config.getfloat('detector','psize')
        self.detector.binned = config.getint('detector','binned')
        self.detector.Nx = config.getint('detector','Nx')
        self.detector.Ny = config.getint('detector','Ny')
        self.detector.gapsize = config.getfloat('detector','gapsize')
        self.detector.gaporientation = config.get('detector','gaporientation')


class Material:
    """
    Material-class of sample-object.
    """
    def __init__(self,parent,**args):
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
            self.massdensity = DICT_massdensity[materialtype]
            self.cH = DICT_atomic_composition[materialtype][0]
            self.cC = DICT_atomic_composition[materialtype][1]
            self.cN = DICT_atomic_composition[materialtype][2]
            self.cO = DICT_atomic_composition[materialtype][3]
            self.cP = DICT_atomic_composition[materialtype][4]
            self.cS = DICT_atomic_composition[materialtype][5]
    
    def get_fX_discrete(self,SF_X):
        """
        get the scattering factor for an element without interpolation.
        Instead it changes the wavelength to closest sampled one.
        """
        e = DICT_physical_constants['e']
        c = DICT_physical_constants['c']
        h = DICT_physical_constants['h']
        ph_energy_eV = c*h/e/self._parent._parent.source.wavelength
        bestmatch = -1
        for line in range(0,pylab.shape(SF_X)[0]):
            if (abs(SF_X[line,0]-ph_energy_eV) < bestmatch) | (bestmatch == -1):
                bestmatch = abs(SF_X[line,0]-ph_energy_eV)
                energy = SF_X[line,0]
                f = SF_X[line,1]
        if abs(energy-ph_energy_eV) > 0.0000001:
           if self._parent._parent._printmode == MODE_PRINT_STD:
               clout = sys.stdout
           if self._parent._parent._printmode == MODE_PRINT_OBJ:
               clout = self._parent._parent.output._tmploglist
           clout.write("Energymismatch = %f eV -> change wavelength to %e m\n" % (energy-ph_energy_eV,c*h/e/energy))
           self._parent._parent.source.wavelength = c*h/e/energy        
        return f

    def get_fX(self,SF_X):
        """
        get the scattering factor for an element through linear interpolation.
        """
        e = DICT_physical_constants['e']
        c = DICT_physical_constants['c']
        h = DICT_physical_constants['h']
        ph_energy_eV = c*h/e/self._parent._parent.source.wavelength
        return pylab.interp(ph_energy_eV,SF_X[:,0],SF_X[:,1])
 
    def get_f_times_n0(self):
        """
        Obtains average atomic scattering factor times average atom number density.
        It uses values relative atomic composition specified in the material-object and the wavelength from the source-object.
        """
        e = DICT_physical_constants['e']
        c = DICT_physical_constants['c']
        h = DICT_physical_constants['h']
        ph_energy_eV = c*h/e/self._parent._parent.source.wavelength
        u = DICT_physical_constants['u']
        elkey_list = []
        c_list = []
        f_list = []
        for key in self.__dict__.keys():
            if key[0] == 'c':
                elkey_list.append(key[1:])
                exec "c_tmp = self." + key
                c_list.append(c_tmp)
        cnorm_array = numpy.array(c_list) / float(sum(c_list))   
        mav = 0
        fav = 0
        for i in range(0,len(elkey_list)):
            # sum up average atom density
            mav = mav + cnorm_array[i]*DICT_atomic_mass[elkey_list[i]]*u
            # sum up average atom factor
            fav += cnorm_array[i]*self.get_fX(DICT_scattering_factors[elkey_list[i]])
        n0 = self.massdensity/mav
        return n0*fav

class DensitymapSample:
    """
    Sample class of input-object that is used to sample by using a densitymap
    """
    def __init__(self,parent,size=None):
        self._parent = parent
        if size != None:
            d = self._parent.source.wavelength*self._parent.detector.distance/self._parent.detector.psize/self._parent.detector._N()
            N = int(round(size/d))
            self.densitymap_d = d
            self.densitymap2d = numpy.zeros((N,N))
            self.densitymap3d = numpy.zeros((N,N,N))
    
    def put_sphere(self,radius,x,y,**materialargs):
        """
        Superpose densitymap of spherical object to 2-dimensional denstiymap (3-dimensional densitymap is not changed).
        Arguments:
        - radius [SI-unit]
        - x,y: sphere positions measured from the center of the densitymap [SI-unit]
        - materialarguments:
          specify a materialtype (example: materialtype='protein',materialtype='virus',materialtype='cell',materialtype='latexball',materialtype='water')
          OR
          relative atomic composition and massdensity (example: cH=2,cO=1,massdensity=1000)

        """
        if self._parent._printmode == MODE_PRINT_STD:
            clout = sys.stdout
        if self._parent._printmode == MODE_PRINT_OBJ:
            clout = self._parent.output._tmploglist

        d = self._parent.source.wavelength*self._parent.detector.distance/self._parent.detector.psize/self._parent.detector._N()
        material_obj = Material(self,**materialargs)
        dm2d_sphere = self._makedm_sphere(radius,material_obj)[0]
        dm2d_balled = self.densitymap2d
        N_balled = len(dm2d_balled)
        N_sphere = len(dm2d_sphere)
        if round(2*(max([x,y])/d+radius/d)) > N_balled:
            N_balled = round(2*(max([x,y])/d+radius/d))
            self.frame(N_balled)                
        for iy in range(0,N_sphere):
            for ix in range(0,N_sphere):
                if self.densitymap2d[int(round(N_balled/2.0+y/d-N_sphere/2.0))+iy,int(round(N_balled/2+x/d-N_sphere/2.0))+ix] != 0 and dm2d_sphere[iy,ix] != 0:
                    clout.write("WARNING: overlap of sphere and given densitymap-object. Further increase density-value of pixel.")
                self.densitymap2d[int(round(N_balled/2.0+y/d-N_sphere/2.0))+iy,int(round(N_balled/2+x/d-N_sphere/2.0))+ix] += dm2d_sphere[iy,ix]

    def _makedm_sphere(self,radius,material_obj):
        """ Creates densitymap of homogeneous sphere """
        if self._parent._printmode == MODE_PRINT_STD:
            clout = sys.stdout
        if self._parent._printmode == MODE_PRINT_OBJ:
            clout = self._parent.output._tmploglist
        f_times_n0 = material_obj.get_f_times_n0()
        d = self._parent.source.wavelength*self._parent.detector.distance/self._parent.detector.psize/self._parent.detector._N()

        r_pix = d*pow(3/(4*numpy.pi),1/3.0)
        N = int(round(2*radius/d))+1
        dm2d = numpy.ones((N,N+1))*N
        dm3d = numpy.ones((N,N,N))

        for iz in range(0,N):
            for iy in range(0,N):
                for ix in range(0,N):
                    delta = numpy.sqrt((ix*d-radius)**2+(iy*d-radius)**2+(iz*d-radius)**2) - radius 
                    if  delta > r_pix:
                        dm2d[iy,ix] -= 1 
                        dm3d[iz,iy,ix] = 0
                    elif delta > -r_pix:
                        dm3d[iz,iy,ix] = 0.5+delta**3/4/r_pix**3
                        dm2d[iy,ix] -= 0.5-delta**3/4/r_pix**3
        dm2d = dm2d*d
        return [f_times_n0*dm2d,f_times_n0*dm3d]

    def _makedm_icosahedron(self,radius,material_obj,eul_ang1,eul_ang2,eul_ang3,speedup_factor):  
        """
        Returns densitymaps [densitymap_2d,densitymap_3d] of homogeneous icosahedron (volume equals volume of a sphere with given radius)
        arguments:
        - radius in m
        - material_object (Intance: Material)
        - 3 Euler-angles for rotation in 3d-space
        - speedup_factor = 1,2,3,... can be set >1 to get a rough densitymap that has a lower resolution as the highest resolution that can be resolved by detector.
        """

        if self._parent._printmode == MODE_PRINT_STD:
            clout = sys.stdout
        if self._parent._printmode == MODE_PRINT_OBJ:
            clout = self._parent.output._tmploglist

        a = radius*(16*numpy.pi/5.0/(3+numpy.sqrt(5)))**(1/3.0)
        Rmax = numpy.sqrt(10.0+2*numpy.sqrt(5))*a/4.0
        Rmin = numpy.sqrt(3)/12*(3.0+numpy.sqrt(5))*a

        f_times_n0 = material_obj.get_f_times_n0()
        d = self._parent.source.wavelength*self._parent.detector.distance/self._parent.detector.psize/self._parent.detector._N()*speedup_factor
 
        nRmax = (Rmax/d)
        nRmin = (Rmin/d) 
        N = int(2*nRmax)

        dm3d = numpy.ones((N,N,N))

        r_pix = d*(3/(4*numpy.pi))**(1/3.0)

        clout.write("... build icosahedron geometry ...\n")
        phi = (1+numpy.sqrt(5))/2.0
        R = 1.0
        x1 = numpy.array([0.0,R,phi])
        x2 = numpy.array([0.0,R,-phi])
        x3 = numpy.array([0.0,-R,phi])
        x4 = numpy.array([0.0,-R,-phi]) 
        x5 = numpy.array([R,phi,0.0])
        x6 = numpy.array([R,-phi,0.0])
        x7 = numpy.array([-R,phi,0.0])
        x8 = numpy.array([-R,-phi,0.0])
        x9 = numpy.array([phi,0.0,R])
        x10 = numpy.array([-phi,0.0,R])
        x11 = numpy.array([phi,0.0,-R])
        x12 = numpy.array([-phi,0.0,-R])
        def cont_element(el,l):
            for i in range(0,len(l)):
                if (el == l[i]).all():
                    return True
            return False
        an = round(numpy.dot(x5,x1))
        def angles_match(y1,y2,y3):
            if round(numpy.dot(y1,y2)) == an and round(numpy.dot(y2,y3)) == an and round(numpy.dot(y3,y1)) == an:
                return True
            else:
                return False
        def rotate_X(v,alpha):
            rotM = numpy.array([[1,0,0],[0,numpy.cos(alpha),-numpy.sin(alpha)],[0,numpy.sin(alpha),numpy.cos(alpha)]])
            return numpy.dot(rotM,v)
        def rotate_Z(v,alpha):
            rotM = numpy.array([[numpy.cos(alpha),-numpy.sin(alpha),0],[numpy.sin(alpha),numpy.cos(alpha),0],[0,0,1]])
            return numpy.dot(rotM,v)
        X = [x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12]
        n_list = []
        for i in range(0,len(X)):
            for j in range(0,len(X)):
                for k in range(0,len(X)):
                    n = (X[i]+X[j]+X[k])/6*a/d
                    if angles_match(X[i],X[j],X[k]) and not cont_element(n,n_list):
                        n_list.append(n)
        clout.write("... build icosahedron in %i x %i x %i grid (%i datapoints)...\n" % (N,N,N,N**3))

        for i in range(0,len(n_list)):
            n_list[i] = rotate_Z(n_list[i],eul_ang1)
            n_list[i] = rotate_X(n_list[i],eul_ang2)
            n_list[i] = rotate_Z(n_list[i],eul_ang3)
        cutpos = []
        for iz in range(0,N):
            for iy in range(0,N):
                for ix in range(0,N):
                    r = (ix-nRmax)**2+(iy-nRmax)**2+(iz-nRmax)**2
                    if r > nRmax**2:
                        dm3d[iz,iy,ix] = 0
                    elif r >= nRmin**2:
                        cutpos.append([iz,iy,ix])
        clout.write("... reduced number of datapoints to %i ...\n" % len(cutpos))
        for m in range(0,len(n_list)):
            n = n_list[m]
            for pos in cutpos:
                r = numpy.array([pos[0]-nRmax,pos[1]-nRmax,pos[2]-nRmax])
                delta = numpy.dot((r-n)*d,1.0*n/nRmin)
                if delta > r_pix:
                    #print "haha"
                    dm3d[pos[0],pos[1],pos[2]] = 0
                elif delta > -r_pix and dm3d[pos[0],pos[1],pos[2]] != 0:
                    #print "huhu"
                    dm3d[pos[0],pos[1],pos[2]] = 0.5+delta**3/4/r_pix**3
                #else: print "ho"
            clout.write("... %i percent done ...\n" % (int(100.0*(m+1)/len(n_list))))
        clout.write("... project icosahedron to plane (%i x %i) ...\n" % (N,N))
        dm2d = self.project(dm3d,d)
        dm3d = dm3d*f_times_n0
        dm2d = dm2d*f_times_n0
        return [dm2d,dm3d,d]

    def project(self,dm3d,dm_d):
        """ Projects 3-dimensional densitymap on 2-dimensional plane """
        N = len(dm3d)
        dm2d = numpy.zeros((N,N))
        for iz in range(0,N):
            for iy in range(0,N):
                for ix in range(0,N):
                    dm2d[iy,ix] = dm2d[iy,ix] + dm3d[iz,iy,ix]*dm_d
        return dm2d

    def frame(self,N_new):
        """ Adds to 2-dimensional densitymap frame (values = zero) in order to enlarge densitymap to N_new x N_new """
        N_old = len(self.densitymap2d)
        densitymap2d_new = numpy.zeros((N_new,N_new))
        for iy in range(0,N_old):
            for ix in range(0,N_old):
                densitymap2d_new[(N_new-N_old)/2+iy,(N_new-N_old)/2+ix] = self.densitymap2d[iy,ix]
        self.densitymap2d = densitymap2d_new

    def pickle_to_file(self,dm,n_dimensions,filename):
        """ Writes 2-dimensional / 3-dimensional densitymap to file """
        file_dm = open(filename,'w')
        if n_dimensions == 2:
            Nx = len(dm[0])
            Ny = len(dm)
            file_dm.write("x y f\n\n")
            for iy in range(0,Ny):
                for ix in range(0,Nx):
                    file_dm.write(" %f\n" % (iy,ix,dm[iy][ix]))
        else:
            N = len(dm)
            for iz in range(0,N):
                for iy in range(0,N):
                    for ix in range(0,N):
                        file_dm.write("%f" % (ix,iy,iz,dm[iz][iy][ix]))      
        file_dm.close()

    def get_area(self,dm,densitymap_d):
        """ Calculates area of object in densitymap """
        area = 0
        for iy in range(0,len(dm)):
            for ix in range(0,len(dm[0])):
                if not dm[iy,ix] == 0:
                    area = area + densitymap_d**2
        return area

    def resize(self,dm2d_old,dm2d_new,N_new):
        """ Resizes given densitymap to new edgelength N_new [datapoints] """
        N_old = len(dm2d_old)
        if N_new < N_old:
            dm2d_new = numpy.zeros(shape=(N_new,N_new))
            n_to_merge = N_old/N_new
            for iy in range(0,N_new):
                for ix in range(0,N_new):
                    for jy in range(0,n_to_merge):
                        for jx in range(0,n_to_merge):
                            dm2d_new[iy,ix] = dm2d_new[iy,ix] + self.dm2d_old[iy*n_to_merge+jy,ix*n_to_merge+jx]/n_to_merge**2
        elif N_new > N_old:
            dm2d_new = numpy.zeros(shape=(N_new,N_new))
            n_to_copy = N_new/N_old
            for iy in range(0,N_new):
                for ix in range(0,N_new):
                    dm2d_new[iy,ix] = self.dm2d_old[iy/n_to_copy,ix/n_to_copy]
        else:
            dm2d_new = dm2d_old

    def plot3d_densitymap(self,densitymap3d=None):
        """
        Creates scatterplot of 3-dimensional densitymap
        If there is no densitymap as argument given function plots self.densitymap3d
        """
        import mpl_toolkits.mplot3d
        if denstiymap3d == None:
            densitymap3d = self.densitymap3d
        N = len(plotmap3d)
        dm3dX = []
        dm3dY = []
        dm3dZ = []
        for iz in range(0,N):
            for iy in range(0,N):
                for ix in range(0,N):
                    if not densitymap3d[iz,iy,ix] == 0:
                        dm3dX.append(ix)
                        dm3dY.append(iy)
                        dm3dZ.append(iz)
        
        fig = pylab.figure()
        #ax = fig.add_subplot(111,projection='3d')
        ax = mpl_toolkits.mplot3d.Axes3D(fig)
        ax.scatter(dm3dX,dm3dY,dm3dZ)
        ax.set_xlabel('X Label')
        ax.set_ylabel('Y Label')
        ax.set_zlabel('Z Label')
        fig.show()

    def plot2d_densitymap(self,plotmap2d=None):
        """
        Creates plot of 2-dimensional densitymap
        If there is no densitymap as argument given function plots self.densitymap2d
        """
        if plotmap2d == None:
            plotmap2d = self.densitymap2d
        fig =  pylab.figure()
        pylab.imshow(plotmap2d,interpolation='nearest')
        fig.show()

class HomogeneoussphereSample:
    """
    Sample class of input-object in order to simulate a homogeneous sphere
    """
    def __init__(self,parent,radius=225E-09,**materialargs):
        self._parent = parent
        self.radius = radius
        self.material = Material(self,**materialargs)

class Detector:
    """ DETECTOR OBJECT
    Initialization of detector: detector_obj = Detector()
    -> creates detector object using default values for pixelsize, gap,..."""
    def __init__(self,parent):
        self.distance = 0.15
        self.psize = 16E-06
        self.binned = 16
        self.Nx = 4096
        self.Ny = 4096
        self.gapsize = 0.0007
        self.gaporientation = 'x'
        self._parent = parent
        
    def _N(self):
        if self.gaporientation == 'x':
            return max([self.Nx,self.Ny+int(self.gapsize/self.psize)])
        if self.gaporientation == 'y':
            return max([int(self.Nx+self.gapsize/self.psize),self.Ny])

    # functions that convert detector-coordinates:
    def _get_q_from_r(self,r):
        return r*2*numpy.pi/self.wavelength/self.distance
    
    def _get_r_from_q(self,q):
        return q/(2*numpy.pi/self._parent.source.wavelength/self.distance)
    
    def _get_i_from_q(self,q):
        """ i is the array-position of value that represents scattered photons to scattering vector q """
        PRES = int(self._N()/self.binned)
        return int(q*self._parent.source.wavelength*PRES*self.distance/(2*numpy.pi*self._N()*self.psize))
   
    def _get_q_from_i(self,i):
        """ i is the array-position of value that represents scattered photons to abs. scattering vector q """
        PRES = int(self._N()/self.binned)
        return i/(self._parent.source.wavelength*PRES*self.distance/(2*numpy.pi*self._N()*self.psize))
    
    def _get_i_from_r(self,r):
        """ i is the array-position of value that represents scattered photons that are scattered to r (spherical coordinates on detector)"""
        PRES = int(self._N()/self.binned)
        return int(r*2*numpy.pi/self._parent.source.wavelength/self.distance*self._parent.source.wavelength*PRES*self.distance/(2*numpy.pi*self._N()*self.psize))


class Source:
    """ SOURCE OBJECT
    Initialization of the source: source_obj = Source()
    -> creates detector object using default values for pixelsize, gap,..."""
    def __init__(self,parent):
        self._parent = parent
        self.wavelength = 5.7E-09 
        self.sizex = 20E-06
        self.sizey = 20E-06
        self.energy = 100E-06

class Output:
    """ Output object of 'spow'
    - contains data generated by 'spow' 
      in order to read date use
      - get_pattern and
      - get_radial_distribution 
    - contains functions for plotting data
      - plot_pattern
      - plot_radial_distribution
    - contains function for saving data
      - save_to_file: saves data to png- or h5-file
      - save_Output_object_to_file: pickles the whole data of the output-object to a file that can be recovered by using load_Output_object_from_file  
    """
    
    def get_pattern(self,scaling="meter"):
        if scaling == "meter":
            return self.intensity_pattern
        elif scaling == "pixel":
            return self.intensity_pattern*self.pixel_size**2
        elif scaling == "nyquist pixel":
            return self.intensity_pattern*self.nyquistpixel_size**2
        else:
            print "ERROR: %s is no valid scaling." % scaling

    def get_radial_distribution(self,scaling="meter",mode="radial average"):
        if mode == "radial average":
            data = self.intensity_radial_average
        elif mode == "radial sum":
            data = self.intensity_radial_sum
        else:
            print "ERROR: %s is no valid mode." % mode
            return
        if scaling == "meter":
            return data
        elif scaling == "pixel":
            return data*self.pixel_size**2
        elif scaling == "nyquist pixel":
            return data*self.nyquistpixel_size**2
        else:
            print "ERROR: %s is no valid scaling." % scaling
            return

    def plot_radial_distribution(self,scaling="pixel and nyquist pixel",mode="all",noise=None):
        """
        Creates 1-dimensional plot(s) showing radial distribution of scattered photons.
        Usage: plot_radial_distribution([scaling],[mode],[noise])
        Arguments:
        - scaling: Specifies spatial scaling.
                   Can be set to 'pixel', 'nyquist pixel', 'pixel and nyquist pixel' or 'meter'.
                   'pixel and nyquist pixel' leads to creation of two plots in one figure using pixel- and Nyquist-pixel-scaling.
        - mode:    Mode specifies whether the radial average or the radial sum will be plotted.
                   Can be set to 'radial average', 'radial sum' or 'all'.
        - noise:   Specifies noise and can be set to 'poisson'.
        """
        if noise == 'poisson':
            def noise(data): return pylab.poisson(data)
        else:
            def noise(data): return data
        def get_arguments(sc):
            if mode == "all":
                legend_args = [('Radial sum', 'Radial average'),'upper right']
                if sc == "pixel":
                    r = numpy.arange(0,len(self.intensity_radial_sum),1)
                elif sc == "nyquist pixel":
                    r = numpy.arange(0,min([self.nyquistpixel_number_x,self.nyquistpixel_number_y])/2,min([self.nyquistpixel_number_x,self.nyquistpixel_number_y])/2/len(self.intensity_radial_sum))
                plot_args = [r,noise(self.get_radial_distribution(sc,'radial sum')),'k',r,noise(self.get_radial_distribution(sc,'radial average')),'k:']
            else:
                if sc == "pixel":
                    r = numpy.arange(0,len(self.intensity_radial_sum),1)
                elif sc == "nyquist pixel":
                    r = numpy.arange(0,min([self.nyquistpixel_number_x,self.nyquistpixel_number_y])/2,min([self.nyquistpixel_number_x,self.nyquistpixel_number_y])/2/len(self.intensity_radial_sum))
                elif sc == "meter":
                    r = numpy.arange(0,min([self.nyquistpixel_number_x,self.nyquistpixel_number_y])/2*self.pixel_size,min([self.nyquistpixel_number_x,self.nyquistpixel_number_y])/2*self.pixel_size/len(self.intensity_radial_sum))
                if mode == "radial sum":
                    legend_args = [('Radial sum'),'upper right']
                    plot_args = [r,noise(self.get_radial_distribution(sc,mode)),'k']
                elif mode == "radial average":
                    legend_args = [('Radial average'),'upper right']
                    plot_args = [r,noise(self.get_radial_distribution(sc,mode)),'k']
            return [plot_args,legend_args]

        if scaling == "pixel and nyquist pixel":
            f1d = pylab.figure(figsize=(10,5))
            f1d.suptitle("\nRadial distribution of scattered photons in detector plane", fontsize=16)
            str_scaling = "binned-pixel"
            f1d_ax_left = f1d.add_axes([0.1, 0.1, 0.35, 0.7],title='Radial scaling:' + str_scaling,xlabel="r [" + str_scaling + "]",ylabel="I(r) [photons/" + str_scaling + "]")
            str_scaling = "Nyquist-pixel"
            f1d_ax_right = f1d.add_axes([0.55, 0.1, 0.35, 0.7],title='Radial scaling:' + str_scaling,xlabel="r [" + str_scaling + "]",ylabel="I(r) [photons/" + str_scaling + "]")
            [plot_args,legend_args] = get_arguments('pixel')
            f1d_ax_left.semilogy(*plot_args)
            f1d_ax_left.legend(*legend_args)
            [plot_args,legend_args] = get_arguments('nyquist pixel')
            f1d_ax_right.semilogy(*plot_args)
            f1d_ax_right.legend(*legend_args)
            f1d.show()
            return
        elif scaling == "pixel":
            str_scaling = "binned pixel"
            r = numpy.arange(0,len(self.intensity_radial_sum),1)
        elif scaling == "nyquist pixel":
            str_scaling == "Nyquist-pixel"
            r = numpy.arange(0,min([self.nyquistpixel_number_x,self.nyquistpixel_number_y])/2,min([self.nyquistpixel_number_x,self.nyquistpixel_number_y])/2/len(self.intensity_radial_sum))
        elif scaling == "meter":
            str_scaling = "meter"
            r = numpy.arange(0,min([self.pixel_number_x,self.pixel_number_y])/2*self.pixel_size,min([self.pixel_number_x,self.pixel_number_y])/2*self.pixel_size/len(self.intensity_radial_sum))
        else:
            print "ERROR: %s is no valid scaling" % scaling
            return
        [plot_args,legend_args] = get_arguments(r,scaling)
        f1d = pylab.figure(figsize=(5,5))
        f1d.suptitle("\nRadial distribution of scattered photons in detector plane", fontsize=16)
        f1d_ax = f1d.add_axes([0.2, 0.1, 0.7, 0.7],title='Radial scaling:' + str_scaling,xlabel="r [" + str_scaling + "]",ylabel="I(r) [photons/" + str_scaling + "]")
        f1d_ax.semilogy(*plot_args)
        f1d_ax.legend(*legend_args)
        f1d.show()
        
    def plot_pattern(self,scaling="pixel and nyquist pixel",noise=None):
        """
        Creates 2-dimensional plot(s) of the distribution of scattered photons.
        Usage: plot_pattern([scaling],[noise])
        Arguments:
        - scaling: Specifies spatial scaling.
                   Can be set to 'pixel', 'nyquist pixel', 'pixel and nyquist pixel' (default) or 'meter'.
                   'pixel and nyquist pixel' leads to creation of two plots in one figure using pixel- and Nyquist-pixel-scaling.
        - noise:   Specifies noise and can be set to 'poisson'.
        """
        if noise == 'poisson':
            def noise(data): return pylab.poisson(data)
        else:
            def noise(data): return data
        if scaling == "pixel and nyquist pixel":
            f2d = pylab.figure(figsize=(10,6))
            # draw intensity plot (N/pixel)
            str_scaling = "binned-pixel"
            max_x = self.pixel_number_x
            max_y = self.pixel_number_y
            f2d.suptitle("\n2-dimensional distribution of scattered photons in detector plane", fontsize=16)
            f2d_ax_left = f2d.add_axes([3/30.0,5/18.0,10/30.0,10/18.0],title='Scaling: ' + str_scaling,xlabel="x [" + str_scaling + "]",ylabel="y [" + str_scaling + "]")
            f2d_axcolor_left = f2d.add_axes([3/30.0,3/18.0,10/30.0,0.5/18.0])
            im_left = f2d_ax_left.matshow(numpy.log10(noise(self.get_pattern('pixel'))),extent=[-max_x/2,max_x/2,-max_y/2,max_y/2])
            cb1 = f2d.colorbar(im_left, cax=f2d_axcolor_left,orientation='horizontal')
            cb1.set_label("log10( I [photons/" + str_scaling + "] )")
            # draw intensity plot (N/Nyquist-pixel)
            str_scaling = "Nyquist-pixel"
            max_x = self.nyquistpixel_number_x
            max_y = self.nyquistpixel_number_y
            f2d_ax_right = f2d.add_axes([17/30.0,5/18.0,10/30.0,10/18.0],title='Scaling: ' + str_scaling,xlabel="x [" + str_scaling + "]",ylabel="y [" + str_scaling + "]")
            f2d_axcolor_right = f2d.add_axes([17/30.0,3/18.0,10/30.0,0.5/18.0])
            im_right = f2d_ax_right.matshow(numpy.log10(noise(self.get_pattern('nyquist pixel'))),extent=[-max_x/2,max_x/2,-max_y/2,max_y/2])
            cb2 = f2d.colorbar(im_right, cax=f2d_axcolor_right,orientation='horizontal')
            cb2.set_label("log10( I [photons/" + str_scaling + "] )")
            f2d.show()
            return
        elif scaling == "meter":
            str_scaling = "meter"
            max_x = self.pixel_number_x*self.pixel_size
            max_y = self.pixel_number_y*self.pixel_size
        elif scaling == "pixel":
            str_scaling = "binned-pixel"
            max_x = self.pixel_number_x
            max_y = self.pixel_number_y
        elif scaling == "nyquist pixel":
            str_scaling = "Nyquist-pixel"
            max_x = self.nyquistpixel_number_x
            max_y = self.nyquistpixel_number_y
        f2d = pylab.figure(figsize=(5,6))
        f2d.suptitle("\n2-dimensional distribution of scattered photons in detector plane", fontsize=16)
        f2d_ax = f2d.add_axes([3/15.0,5/18.0,10/15.0,10/18.0],title='Scaling: ' + str_scaling,xlabel="x [" + str_scaling + "]",ylabel="y [" + str_scaling + "]")
        f2d_axcolor = f2d.add_axes([3/15.0,3/18.0,10/15.0,0.5/18.0])
        im = f2d_ax.matshow(numpy.log10(noise(self.get_pattern(scaling))),extent=[-max_x/2,max_x/2,-max_y/2,max_y/2])
        cb = f2d.colorbar(im, cax=f2d_axcolor,orientation='horizontal')
        cb.set_label("log10( I [photons/" + str_scaling + "] )")
        f2d.show()
   
    def save_pattern_to_file(self,filename,scaling="pixel",*arguments):
        """
        Saves dataset to file of specified format.
        Usage: fo_file(filename,[scaling],[colorscale])
        Arguments:
        - filename: The file-format is specified using one of the following file-endings:
                    - '.h5'
                    - '.png'
        - scaling:  Specifies spatial scaling.
                    Can be set to 'pixel' (default), 'nyquist pixel' or 'meter'.
        - colorscale (only for png-files):
                    - Jet
                    - Gray (default)
                    - Log (can be combined with the others)
        """
        import spimage,h5py
        pattern = self.get_pattern(scaling)
        if filename[-3:]=='.h5':
            color = 0
        elif filename[-3:]=='.png':
            color = 16
            for flag in arguments:
                if flag == 'Jet':
                    color = 16
                elif flag == 'Gray':
                    color = 1
                elif flag == 'Log':
                    color += 128
                else:
                    print "unknown flag %s" % flag
                    return
        else:
            print "ERROR: %s is not a valid fileformat for this function." % filename[-3:]
            return
        tmp_data = spimage.sp_image_alloc(len(pattern[0]),len(pattern),color)
        tmp_data.image[:,:] = pattern[:,:]
        spimage.sp_image_write(tmp_data,filename,0)
        spimage.sp_image_free(tmp_data)

    def save_Output_object_to_file(self,filename):
        picklefile = open(filename,'w')
        keys = self.__dict__.keys()
        topickle_list = [keys]
        for key in keys:
            exec "tmp = self." + key
            topickle_list.append(tmp)
        pickle.dump(topickle_list,picklefile)
        picklefile.close()
    
def load_Output_object_from_file(filename):
    unpicklefile = open(filename,'r')
    unpickled_list = pickle.load(unpicklefile)
    keys = unpickled_list[0]
    output_obj = Output()
    for i in range(0,len(keys)):
        exec "output_obj." + keys[i] + " = unpickled_list[i+1]"
    return output_obj        

# Class for commandline-output
#-----------------------------
class _WritableObject:
    """ Class that can be assigned to stdout in order to switch off output to commandline and instead to write it to self.content"""
    def __init__(self):
        self.content = []
    def write(self, string):
        self.content.append(string)


# MAIN FUNCTION
#--------------
def spow(input_obj=False):
    """ MAIN FUNCTION of 'spowpy.py': 'spow' calculates distribution of scattered photons by defined sample.
    Usage: output_obj = spow(input_obj)
    input_obj: configuration object of class 'Input' that sets simulation and program parameters
    output_obj: contains data generated by 'spow' and plotting-functions"""

    # USAGE
    # -----

    # if given input objects incorrect print usage information
    if not isinstance(input_obj,Input):
        print "... ERROR: WRONG INPUT ...\n" 
        print "Usage: spow(input_obj)\n"
        print "   input_obj: configuration object that sets simulation parameters (instance: Input)"
        return
    
    # INITIALIZATION
    # --------------
    
    # set modes
    plotmode = input_obj._plotmode
    printmode = input_obj._printmode
    output = Output()
    if printmode == MODE_PRINT_STD:
        clout = sys.stdout
    if printmode == MODE_PRINT_OBJ:
        clout = input_obj.output._tmploglist
 
    # read physical constants from dictionary (SI units):
    e = DICT_physical_constants['e']
    c = DICT_physical_constants['c']
    h = DICT_physical_constants['h']
    re = DICT_physical_constants['re']
    barn = DICT_physical_constants['barn']

    # read / calculate experimental paramters   
    ph_energy_eV = c*h/e/input_obj.source.wavelength
    if ph_energy_eV > F_MAX_ENERGY_EV or ph_energy_eV < F_MIN_ENERGY_EV:
        print ph_energy_eV
        print "ERROR: Photons at the given energy can't be simulated. Appropriate atomic scattering factors are missing."
        print "       Minimal photon energy = %f eV" % F_MIN_ENERGY_EV
        print "       Maximal photon energy = %f eV" % F_MAX_ENERGY_EV
        return

    # abbreviations
    so_wavelength = input_obj.source.wavelength 
    ph_energy_eV = c*h/e/so_wavelength
    so_sizex = input_obj.source.sizex
    so_sizey = input_obj.source.sizey
    so_energy = input_obj.source.energy
    ph_energy = c*h/so_wavelength
    so_photons = so_energy/ph_energy
    so_area = input_obj.source.sizex*input_obj.source.sizey
    I0 = so_photons/so_area
    de_distance = input_obj.detector.distance
    de_psize = input_obj.detector.psize * input_obj.detector.binned
    de_Nx_binned = input_obj.detector.Nx / input_obj.detector.binned
    de_Ny_binned = input_obj.detector.Ny / input_obj.detector.binned
    de_gapsize = input_obj.detector.gapsize
    de_gaporientation = input_obj.detector.gaporientation
    if de_gaporientation == 'x':
        de_Ny_binned = de_Ny_binned + int(de_gapsize/de_psize)
    elif de_gaporientation == 'y':
        de_Nx_binned = de_Nx_binned + int(de_gapsize/de_psize)
    de_N_binned_long = max([de_Nx_binned,de_Ny_binned])
    de_N_binned_short = min([de_Nx_binned,de_Ny_binned])
    de_pdOmega = 1.0/de_distance**2    
    dm_d_detector = so_wavelength*de_distance/(de_psize*de_N_binned_long)
    dm_dA_detector = dm_d_detector**2

    # DEFINE GENERAL FUNCTIONS
    # ------------------------
     
    # calculate area [pixel] of object from 2-dimensional densitymap
    def get_supparea(dm2d,dm_d):
        res = 0
        for iy in range(0,len(dm2d[0])):
            for ix in range(0,dm2d):
                if not dm2d[iy,ix] == 0:
                    res = res + 1
        return res*dm_d**2
    
    # calculates 2d distribution of scattered photons using radial function of scattering amplitude 
    def N_2dpattern_radial(func):
        N = numpy.zeros((de_N_binned_long,de_N_binned_long))
        for iy in range(0,de_N_binned_long):
            for ix in range(0,de_N_binned_long):
                N[iy,ix] = func(numpy.sqrt(input_obj.detector._get_q_from_i(iy-de_N_binned_long/2)**2+input_obj.detector._get_q_from_i(ix-de_N_binned_long/2)**2))
        return N
    
    # calculate radial sum and radial average from given 2-dimensional pattern
    def N_radial(pattern):
        Nsum = numpy.zeros(de_N_binned_short/2)
        Nav = numpy.zeros(de_N_binned_short/2)
        for iy in range(0,de_Ny_binned):
            for ix in range(0,de_Nx_binned):
                ir = int(numpy.sqrt((de_Ny_binned/2-iy)**2+(de_Nx_binned/2-ix)**2))
                if ir < de_N_binned_short/2:
                    Nsum[ir] = Nsum[ir] + pattern[iy][ix]
                    Nav[ir] = Nav[ir] + 1
        return [Nsum,Nsum/Nav]

    # calculate absolute number of scattered photons from given pattern
    def abs_scat_photons(Npattern):
        res = 0
        for iy in range(0,len(Npattern)):
            for ix in range(0,len(Npattern[0])):
                res = res + Npattern[iy,ix]*de_N_binned_long**2/de_N_binned_long**2
        return res
    
    # delete gap between detector halves
    def N_delete_gap(Npattern):
        if de_gapsize > 0:
            if de_gaporientation == 'x':
                ix_list = range(0,de_N_binned_long)
                iy_list = range(de_N_binned_long/2-input_obj.detector._get_i_from_r(de_gapsize/2.0),de_N_binned_long/2+input_obj.detector._get_i_from_r(de_gapsize/2.0))
            elif de_gaporientation == 'y':
                ix_list = range(de_N_binned_long/2-input_obj.detector._get_i_from_r(de_gapsize/2.0),de_N_binned_long/2+input_obj.detector._get_i_from_r(de_gapsize/2.0))
                iy_list = range(0,de_N_binned_long)
            for iy in iy_list:
                for ix in ix_list:
                    Npattern[iy,ix] = 0

    # cut border of 2d-pattern to real detector size
    def N_cut_to_real_size(Npattern):
        if de_N_binned_long != de_Nx_binned or de_N_binned_long != de_Ny_binned:
            Npattern_new = numpy.zeros((de_Ny_binned,de_Nx_binned))
            dx = de_N_binned_long-de_Nx_binned
            dy = de_N_binned_long-de_Ny_binned
            for iy in range(0,de_Ny_binned):
                for ix in range(0,de_Nx_binned):
                    Npattern_new[iy][ix] = Npattern[dy/2+iy][dx/2+ix]
            return Npattern_new
        else:
            return Npattern
 

    if isinstance(input_obj.sample,HomogeneoussphereSample):
        # SAMPLEMODE: HOMOGENEOUSSPHERE
        
        # GET SPHERICAL SAMPLE PARAMETERS
        # -------------------------------
        clout.write("... simulate sphere ...\n") 
   
        # read sample radius
        sa_radius = input_obj.sample.radius

        # calculate average atom density and average atomic scattering factor
        f_times_n0 = input_obj.sample.material.get_f_times_n0()
     
        # DEFINE SPECIAL FUNCTIONS
        # ------------------------

        # Scattering amplitude by homogeneous sphere
        def F_homsphere(q):
            if q==0:
                return 0
            else:
                return 4*numpy.pi*f_times_n0*(numpy.sin(q*sa_radius)-q*sa_radius*numpy.cos(q*sa_radius))/pow(q,3)

        # Scattered photons per pixel
        def N_homsphere_pix(q):
            return I0*de_pdOmega*pow(re,2)*pow(F_homsphere(q),2)

        # BUILD PATTERN
        # -------------
        clout.write("... build %s x %s pattern ...\n" % (de_N_binned_long,de_N_binned_long))
 
        # Create pattern de_N_binned_long x de_N_binned_long
        Npattern = N_2dpattern_radial(N_homsphere_pix)

        # Calculate Nyquist pixelsize considering area of sample
        dq_nyquist = numpy.pi/sa_radius
        dr_nyquist = input_obj.detector._get_r_from_q(dq_nyquist)
        max_x_Nypix = de_Nx_binned*de_psize/dr_nyquist 
        max_y_Nypix = de_Ny_binned*de_psize/dr_nyquist 

   
    elif isinstance(input_obj.sample,DensitymapSample):
        # SAMPLEMODE: DENSITYMAP
        
        # READ DENSITYMAP
        # ---------------
        clout.write("... simulate density-map-object ...\n")

        # check given densitymap
        dm_d_given = input_obj.sample.densitymap_d
        dm_dA_given = dm_d_given**2
        try:
            # 3d
            dm3d = input_obj.sample.densitymap3d
            clout.write("... 3-dimensional density map given ...\n")
            dm_dim = 3
            # apply projection approximation
            clout.write("... project 3-dimensional density map to plane ...\n")
            input_obj.sample.project(dm3d,dm2d,dm_d_given)
        except:
            # 2d
            clout.write("... 2-dimensional density map given ...\n")
            dm_dim = 2
            dm2d = input_obj.sample.densitymap2d
      
        # resize if necessary
        N_new = int(round(dm_d_given/dm_d_detector*len(dm2d)))
        input_obj.sample.resize(dm2d,dm2d,N_new)
        
        #input_obj.sample.plot2d_densitymap(dm2d)

        # DEFINE SPECIAL FUNCTIONS
        # ------------------------
 
        # Scattered photons per pixel
        def N_dm_2dpattern_pix():
            return I0*re**2*(dm_dA_detector*abs(numpy.fft.fftn(dm2d,(de_N_binned_long,de_N_binned_long))))**2*de_pdOmega

        # BUILD PATTERN
        # -------------
        clout.write("... build %s x %s pattern ...\n" % (de_N_binned_long,de_N_binned_long))    


        # Create pattern de_N_binned_long x de_N_binned_long
        Npattern = pylab.fftshift(N_dm_2dpattern_pix())

        # Calculate Nyquist pixelsize considering area of sample
        sa_area = input_obj.sample.get_area(dm2d,dm_d_detector)
        assumed_sa_radius = numpy.sqrt(sa_area/numpy.pi)
        dq_nyquist = numpy.pi/numpy.sqrt(sa_area/numpy.pi)
        dr_nyquist = input_obj.detector._get_r_from_q(dq_nyquist)
        max_x_Nypix = de_Nx_binned*de_psize/dr_nyquist
        max_y_Nypix = de_Ny_binned*de_psize/dr_nyquist
    
    else:
        print "ERROR: No valid sample-object given."
        return

    # Delete gaps between detector halves and cut to real size
    N_delete_gap(Npattern)
    Npattern = N_cut_to_real_size(Npattern)

    # Determine radial sums of scattered photons
    [Nradsum,Nradav] = N_radial(Npattern)


    # COPY OUTPUT TO OUTPUT-OBJECT
    # ----------------------------
   
    clout.write("... write output-object ...\n")

    output.input_object = input_obj

    output.wavelength = so_wavelength
    output.scattered_photons_detector =  abs_scat_photons(Npattern)

    output.intensity_pattern = Npattern
    output.intensity_radial_sum = Nradsum
    output.intensity_radial_average = Nradav

    output.pixel_size = de_psize
    output.nyquistpixel_size = dr_nyquist
    output.pixel_number_x = de_Nx_binned
    output.pixel_number_y = de_Ny_binned
    output.nyquistpixel_number_x = max_x_Nypix
    output.nyquistpixel_number_y = max_y_Nypix


    # PLOTTING
    # --------
 
    if  (plotmode & MODE_PLOT_1D) != 0:
        # plot radial distribution of intensity pattern
        clout.write("... generate 1d-plot ...\n")
        output.plot_radial_distribution()
                
    if (plotmode & MODE_PLOT_2D) != 0:
        # plot 2-dimensional distribution of intensity pattern
        clout.write("... generate 2d-plot ...\n")
        output.plot_pattern()

    # COMMAND LINE OUTPUT
    #--------------------
    if printmode == MODE_PRINT_OBJ:
        # write commandline output to loglist   
        output.loglist = input_obj.output._tmploglist
        input_obj.output._tmploglist = _WritableObject()

    return output

def transform_density(density,density_d,oversampling,limit_resolution = 1):
    """Fourier transforms a density like the one found in
    Sample.densitymap3d."""
    if limit_resolution != 1:
        x = pylab.arange(0,pylab.shape(density)[0],limit_resolution)
        y = pylab.arange(0,pylab.shape(density)[1],limit_resolution)
        z = pylab.arange(0,pylab.shape(density)[2],limit_resolution)
        density = density[x][:,y][:,:,z]
    pattern = ((density_d*limit_resolution)**3*abs(pylab.fftshift(pylab.fftn(density,pylab.array(pylab.shape(density))*oversampling))))**2
    pixel_size = 1.0/((density_d*limit_resolution)*pylab.shape(pattern)[0])

    return [pattern, pixel_size]

def sample_3d_density(density,density_d,source,detector,rotation):
    """We assume that the pixel size in the density is
    the same as the maximum resolution allowed by the
    detector setup"""
    #create 3 arrays with fourier coordinates
    x = pylab.arange(-detector.Nx/2,detector.Nx/2,detector.binned)*detector.psize
    y = pylab.arange(-detector.Ny/2,detector.Ny/2,detector.binned)*detector.psize
    x_angle = pylab.arctan2(x,detector.distance)
    y_angle = pylab.arctan2(y,detector.distance)
    x_fourier = pylab.sin(x_angle)/source.wavelength
    y_fourier = pylab.sin(y_angle)/source.wavelength
    X,Y = pylab.meshgrid(x_fourier,y_fourier)
    Z = (1.0 - pylab.sqrt(1.0 - (X**2+Y**2)*source.wavelength**2))/source.wavelength

    #this should be updated later (transform from fourier coordinates to pixels)
    # X_pixels = X*detector.distance/detector.Nx/detector.psize*wavelength*2
    # Y_pixels = Y*detector.distance/detector.Ny/detector.psize*wavelength*2
    # Z_pixels = Z*detector.distance/detector.Nx/detector.psize*wavelength*2 #bad
    X_pixels = X/density_d
    Y_pixels = Y/density_d
    Z_pixels = Z/density_d

    coord_matrix = pylab.matrix([X_pixels.flatten(),Y_pixels.flatten(),Z_pixels.flatten()])

    #rotate coordinates
    c1 = pylab.cos(rotation[0]); c2 = pylab.cos(rotation[1]); c3 = pylab.cos(rotation[2])
    s1 = pylab.sin(rotation[0]); s2 = pylab.sin(rotation[1]); s3 = pylab.sin(rotation[2])
    rot_mat = pylab.matrix([[c1*c3-c2*s1*s3, -c1*s3-c3*c2*s1, s2*s1],
                            [c2*c1*s3+c3*s1, c1*c2*c3-s1*s3, -c1*s2],
                            [s3*s2, c3*s2, c2]])
    pixels_rot = rot_mat*coord_matrix
    pixels_rot += pylab.outer(pylab.array(pylab.shape(density))/2,
                              pylab.ones(pylab.shape(pixels_rot)[1]))
    #pixels_rot += pylab.outer([detector.Nx/detector.binned/2,detector.Ny/detector.binned/2,detector.Nx/detector.binned/2],pylab.ones(pylab.shape(pixels_rot)[1]))
    #sample

    values = pylab.zeros(pylab.shape(pixels_rot)[1])
    cases = [pylab.floor,lambda x: 1.0 + pylab.floor(x)]
    try:
        for i in range(len(values)):
            if min(pixels_rot[:,i]) >= 0.0 and max(pixels_rot[:,i]) <= pylab.shape(density)[0]-1: #bad
                for x_func in cases:
                    for y_func in cases:
                        for z_func in cases:
                            values[i] +=\
                                (1.0-abs(x_func(pixels_rot[0,i]) - pixels_rot[0,i]))*\
                                (1.0-abs(y_func(pixels_rot[1,i]) - pixels_rot[1,i]))*\
                                (1.0-abs(z_func(pixels_rot[2,i]) - pixels_rot[2,i]))*\
                                density[int(x_func(pixels_rot[0,i])),
                                        int(y_func(pixels_rot[1,i])),
                                        int(z_func(pixels_rot[2,i]))]
    except:
        print "%d %d %d" % (pixels_rot[0,i],pixels_rot[1,i],pixels_rot[2,i])

    print pylab.shape(values)
    print (detector.Nx,detector.Ny)
    
    #rescale depending on solid angle of pixel
    total_angle = pylab.arctan2(pylab.sqrt(x**2+y[:,pylab.newaxis]**2),detector.distance)
    #solid_angle = (detector.psize*detector.binned/detector.distance*pylab.cos(total_angle))**2
    solid_angle = (detector.psize*detector.binned)**2/detector.distance**2*pylab.cos(total_angle)**3
    
    print solid_angle

    values = values.reshape((detector.Nx/detector.binned,detector.Ny/detector.binned))*solid_angle

    #beam intensity and electron cross section
    I0 = source.energy/(DICT_physical_constants['c']*DICT_physical_constants['h']/source.wavelength)/source.sizex/source.sizey
    print "I0 = ",I0
    print "re = ",DICT_physical_constants['re']
    values = I0*DICT_physical_constants['re']**2*values
    #something else?

    #print sys.argv[3]

    return values

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print "... ERROR: WRONG INPUT ...\n"
        print "Usage: python spowpy.py [configuration-file] [plotmode]\n"
        print "   configuration-file: please edit configuration-file to specify detector and sample\n"
        print "   plotmode:   - 'none' -> no plots"
        print "               - '1d'   -> 1d-plots"
        print "               - '2d'   -> 2d-plots"
    else:
        input_obj = Input(sys.argv[1])
        if sys.argv[2] == '1d': input_obj.set_plotmode('1d')
        elif sys.argv[2] == '2d': input_obj.set_plotmode('2d')
        elif sys.argv[2] == 'none': input_obj.set_plotmode('none')
        else: "ERROR: %s is not a valid plotmode." % sys.argv[2]
        spow(input_obj)
        if input_obj._plotmode != MODE_PLOT_NONE:
            pylab.show()
