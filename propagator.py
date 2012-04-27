# ----------------------------------------------------------------------------------------------------- 
# PROPAGATOR: Scattering experiment simulator for spheres and customized object maps
# Please type 'help propagator()' for further information.
# -----------------------------------------------------------------------------------------------------
# Author:  Max Hantke - maxhantke@gmail.com
# -----------------------------------------------------------------------------------------------------

# TO DO:
# - Warn user if Fraunhofer approximation or Born approximation is not holding.
# - Polarization factor needed?

import pylab, sys, ConfigParser, numpy, types, pickle, time, math
from matplotlib import rc
import matplotlib.pyplot as mpy
rc('text', usetex=True)
rc('font', family='serif')

import config
config.init_configuration()

import xcorepropagation,imgutils,tools

from source import *
from sample import *
from detector import *

def propagator(input_obj=False):
    """
    Function calculates diffraction under defined conditions specified in the input object.
    =======================================================================================
    Arguments:
    input_obj: Input object of class \'Input\'.

    """
    
    # Check input
    if not isinstance(input_obj,Input):
        print "ERROR: Illegal input. Argument has to be of instance Input." 
        return
    
    t_start = time.time()

    # Auxiliary variables
    wavelength = input_obj.source.photon.get_wavelength()
    I_0 = input_obj.source.pulse_energy / input_obj.source.photon._energy / input_obj.source.get_area() 
    Omega_p = (input_obj.detector.pixel_size*input_obj.detector.binning)**2 / input_obj.detector.distance**2

    if isinstance(input_obj.sample,SampleSphere):    
        # scattering amplitude from homogeneous sphere: F = sqrt(I_0 Omega_p) 2pi/wavelength^2 [ 4/3 pi R^3  3 { sin(qR) - qR cos(qR) } / (qR)^3 ] dn_real
        dn_real = (1-input_obj.sample.material.get_n()).real
        R = input_obj.sample.radius
        q = input_obj.generate_absqmap()
        F = pylab.zeros_like(q)
        F[q!=0.0] = (pylab.sqrt(I_0*Omega_p)*2*pylab.pi/wavelength**2*4/3.0*pylab.pi*R**3*3*(pylab.sin(q*R)-q*R*pylab.cos(q*R))/(q*R)**3*dn_real)[q!=0.0]
        try: F[q==0] = pylab.sqrt(I_0*Omega_p)*2*pylab.pi/wavelength**2*4/3.0*pylab.pi*R**3*dn_real
        except: pass

    if isinstance(input_obj.sample,SampleMap):    
        # scattering amplitude from dn-map: F = sqrt(I_0 Omega_p) 2pi/wavelength^2 [ DFT{dn_perp} ] dA
        q = input_obj.generate_qmap()
        # requested q interval
        QI = 2*pylab.pi/(input_obj.sample.dX*input_obj.propagation.rs_oversampling)
        # q map scaled for nfft (standard interval -0.5 to 0.5)
        q_scaled = q/QI
        config.OUT.write("Propagate pattern of %i x %i pixels.\n" % (q.shape[1],q.shape[0]))
        F = pylab.sqrt(I_0*Omega_p)*2*pylab.pi/wavelength**2 \
            * xcorepropagation.nfftXCore(input_obj.sample.map3d,
                                         q_scaled,
                                         self.input_object.propagation.N_processes) \
            * input_obj.sample.dX**3
        config.OUT.write("Got pattern of %i x %i pixels.\n" % (F.shape[1],F.shape[0]))
        #F = imgutils.crop(F,pylab.array([input_obj.detector.mask.shape[0],input_obj.detector.mask.shape[1]]))

    t_stop = time.time()
    config.OUT.write("Time = %f sec\n" % (t_stop-t_start))
    config.OUT.write("Propagation finished.\n")
    return Output(F,input_obj)

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
        self.source = Source(parent=self)
        self.sample = SampleSphere(parent=self)
        self.detector = Detector(parent=self)
        self.propagation = Propagation(parent=self)
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

    def load_sample_map(self,filename):
        """
        Creates empty densitymap. Size in accordance to given detector geometry.
        Densitymap resolution is set according to the detector geometry.
        """
        self.sample = SampleMap(parent=self)
        self.sample.load_map3d(filename)

    def set_sample_icosahedral_virus_map(self,radius=None,eul_ang1=0.0,eul_ang2=0.0,eul_ang3=0.0):
        """
        Creates refractive index map of homogeneously filled icosahedron.
        - material: 'virus' (for details investigate generated material object)
        - volume corresponds to sphere of given radius (without given radius set to radius specified in current sample object)
        - 3 euler angles define orientation of the icosahedron in the mesh. (0,0,0) means that the axis of the beam coincides with one of the 2-fold axes
        """
        if not radius:
            radius = self.sample.radius
        self.sample = SampleMap(parent=self)
        self.sample.put_icosahedral_virus(radius,eul_ang1,eul_ang2,eul_ang3)
        self.sample.radius = radius

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
          set_sample_homogeneous_sphere(radius,massdensity=1000,cH=2,cO=1)
        - setting atomic composition values and massdensity according to given materialtype:
          set_sample_homogeneous_sphere(radius,materialtype='protein')
        available materialtypes: 'protein', 'virus', 'cell', 'latexball', 'water', 'Au'
        """ 
        self.sample = SampleSphere(parent=self,**kwargs)

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
        mat = C.get('sample','material')
        args = []
        if mat == 'custom':
            cX_list = C.items('sample')
            for cX_pair in cX_list:
                if cX_pair[0][0] == 'c':
                    el = cX_pair[0]
                    el = el[1:].capitalize()
                    val = float(cX_pair[1])
                    args.append(("'c%s'" % el,val))
            args.append(('massdensity',C.getfloat('sample','massdensity')))
        else:
            keys = ['cH','cN','cO','cP','cS']
            for i in range(0,len(keys)):
                args.append((keys[i],config.DICT_atomic_composition[mat][i]))
            args.append(('massdensity',config.DICT_massdensity[mat]))
        args= dict(args)
        args['radius']=C.getfloat('sample','radius')
        self.set_sample_homogeneous_sphere(**args)
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
        else:
            M = pylab.imread(C.get('detector','mask'))
            M = M[:,:,0]
            M[abs(M)>0] = 1
            M[M==0] = 0
            args['mask'] = M.copy()
        self.detector = Detector(**args)
        self.propagation.rs_oversampling = C.getfloat('propagation','rs_oversampling')
        self.propagation.N_processes = C.getint('propagation','N_processes')

    def get_real_space_resolution_element(self):
        """
        Function returns real space sampling.
        =====================================
        Formula: dX = 2 * pi / max([qxmax,qymax])

        """
        dX = 2 * pylab.pi / max([self.get_absqx_max(),self.get_absqy_max()])
        return dX

    def get_max_achievable_crystallographic_resolution(self):
        """
        Function returns maximal horizontal and vertical crystallographoc resolution achievable in the current setup.
        =============================================================================================================
        Formula: Rx_c = dx * 2 = 4 * pi / qxmax ; Ry_c = dy * 2 = 4 * pi / qymax

        """
        dx = 4 * pylab.pi / self.get_absqx_max()
        dy = 4 * pylab.pi / self.get_absqy_max()
        return [dx,dy]

    def generate_absqmap(self):
        wavelength = self.source.photon.get_wavelength()
        X,Y = pylab.meshgrid(pylab.arange(0,self.detector.mask.shape[1]),pylab.arange(0,self.detector.mask.shape[0]))
        X -= self.detector.cx
        Y -= self.detector.cy
        theta = pylab.arctan2(1.0*Y,1.0*X)
        phi = pylab.arctan2(pylab.sqrt(X**2+Y**2)*self.detector.pixel_size,self.detector.distance)
        R_Ewald = 2*pylab.pi/wavelength
        qmap = R_Ewald * pylab.sqrt((pylab.cos(phi)-1.0)**2 + (pylab.sin(phi) * pylab.sin(theta))**2 + (pylab.sin(phi) * pylab.cos(theta))**2)
        return qmap

    def generate_qmap(self):
        wavelength = self.source.photon.get_wavelength()
        qmap = pylab.zeros(shape=(self.detector.mask.shape[0],
                                  self.detector.mask.shape[1],
                                  3))
        X,Y = pylab.meshgrid(pylab.arange(0,self.detector.mask.shape[1]),pylab.arange(0,self.detector.mask.shape[0]))
        X -= self.detector.cx
        Y -= self.detector.cy
        theta = pylab.arctan2(Y,X)
        phi = pylab.arctan2(pylab.sqrt(X**2+Y**2)*self.detector.pixel_size,self.detector.distance)
        R_Ewald = 2*pylab.pi/wavelength
        qmap[:,:,0] = R_Ewald * (pylab.cos(phi)-1.0)
        qmap[:,:,1] = R_Ewald * pylab.sin(phi) * pylab.sin(theta)
        qmap[:,:,2] = R_Ewald * pylab.sin(phi) * pylab.cos(theta)
        try: E0 = self.sample.euler_angle_0
        except: E0 = 0.0
        try: E1 = self.sample.euler_angle_1
        except: E1 = 0.0
        try: E2 = self.sample.euler_angle_2
        except: E2 = 0.0
        if E0 != 0.0 or E1 != 0.0 or E2 != 0.0:
            M = pylab.array([[pylab.cos(E1)*pylab.cos(E2),
                              -pylab.cos(E0)*pylab.sin(E2)+pylab.sin(E0)*pylab.sin(E1)*pylab.cos(E2),
                              pylab.sin(E0)*pylab.sin(E2)+pylab.cos(E0)*pylab.sin(E1)*pylab.cos(E2)],
                             [pylab.cos(E1)*pylab.sin(E2),
                              pylab.cos(E0)*pylab.cos(E2)+pylab.sin(E0)*pylab.sin(E1)*pylab.sin(E2),
                              -pylab.sin(E0)*pylab.cos(E2)+pylab.cos(E0)*pylab.sin(E1)*pylab.sin(E2)],
                             [-pylab.sin(E1),
                               pylab.sin(E0)*pylab.cos(E1),
                               pylab.cos(E0)*pylab.cos(E1)]])
            for iy in pylab.arange(0,qmap.shape[0]):
                for ix in pylab.arange(0,qmap.shape[0]):
                    qmap[iy,ix,:] = pylab.dot(M,qmap[iy,ix,:])
        return qmap
        
    def get_absqx_max(self):
        wavelength = self.source.photon.get_wavelength()
        x_max = max([self.detector.cx,self.detector.Nx-self.detector.cx])
        R_Ewald = 2*pylab.pi/wavelength
        phi = pylab.arctan2(x_max,self.detector.distance)
        return R_Ewald * pylab.sin(phi)

    def get_absqy_max(self):
        wavelength = self.source.photon.get_wavelength()
        y_max = max([self.detector.cy,self.detector.Ny-self.detector.cy])
        R_Ewald = 2*pylab.pi/wavelength
        phi = pylab.arctan2(y_max,self.detector.distance)
        return R_Ewald * pylab.sin(phi)

class Output:
    """
    OUTPUT of propagator provides user with results and functions for plotting.
    """
    def __init__(self,amplitudes,input_object):
        self.amplitudes = amplitudes.copy()
        self.input_object = input_object 
    
    def get_intensity_pattern(self):
        """
        Returns 2-dimensional array with intensity values in photons per pixel (binned).
        """
        return abs(self.amplitudes)**2

    def get_intensity_radial_average(self):
        """
        Returns 1-dimensional array with intensity average in photons per pixel (binned). x-coordinate sampling is pixel (binned). 
        """
        I = self.get_intensity_pattern()
        return imgutils.radial_pixel_average(I)


    def get_intensity_radial_sum(self):
        """
        Returns 1-dimensional array with intensity average in photons per pixel (binned). x-coordinate sampling is pixel (binned). 
        """
        I = self.get_intensity_pattern()
        return imgutils.radial_pixel_sum(I)

            
    def plot_radial_distribution(self,scaling="binned pixel and nyquist pixel",mode="all",noise=None):
        """
        Creates 1-dimensional plot(s) showing radial distribution of scattered photons.
        Usage: plot_radial_distribution([scaling],[mode],[noise])
        Arguments:
        - scaling: Specifies spatial scaling.
                   Can be set to 'binned pixel', 'nyquist pixel', 'binned pixel and nyquist pixel' or 'meter'.
                   'binned pixel and nyquist pixel' leads to creation of two plots in one figure using pixel- and Nyquist-pixel-scaling.
        - mode:    Mode specifies whether the radial average or the radial sum will be plotted.
                   Can be set to 'radial average', 'radial sum' or 'all'.
        - noise:   Specifies noise and can be set to 'poisson'.
        """
        Ir_avg = self.get_intensity_radial_average()
        Ir_sum = self.get_intensity_radial_sum()
        if noise == 'poisson':
            def noise(data): return pylab.poisson(data)
        else:
            def noise(data): return data
        def get_arguments(sc):
            if mode == "all":
                legend_args = [('Radial sum', 'Radial average'),'upper right']
                if sc == "binned pixel":
                    r = numpy.arange(0,len(Ir_sum),1)
                elif sc == "nyquist pixel":
                    r = numpy.arange(0,min([self.nyquistpixel_number_x,self.nyquistpixel_number_y])/2,min([self.nyquistpixel_number_x,self.nyquistpixel_number_y])/2/len(Ir_sum))
                plot_args = [r,noise(Ir_sum),'k',r,noise(Ir_avg),'k:']
            else:
                if sc == "binned pixel":
                    r = numpy.arange(0,len(Ir_sum),1)
                elif sc == "nyquist pixel":
                    r = numpy.arange(0,min([self.nyquistpixel_number_x,self.nyquistpixel_number_y])/2,min([self.nyquistpixel_number_x,self.nyquistpixel_number_y])/2/len(Ir_sum))
                elif sc == "meter":
                    r = numpy.arange(0,min([self.nyquistpixel_number_x,self.nyquistpixel_number_y])/2*self.pixel_size,min([self.nyquistpixel_number_x,self.nyquistpixel_number_y])/2*self.pixel_size/len(Ir_sum))
                if mode == "radial sum":
                    legend_args = [('Radial sum'),'upper right']
                    plot_args = [r,noise(Ir_sum),'k']
                elif mode == "radial average":
                    legend_args = [('Radial average'),'upper right']
                    plot_args = [r,noise(Ir_avg),'k']
            return [plot_args,legend_args]

        if scaling == "binned pixel and nyquist pixel":
            f1d = pylab.figure(figsize=(10,5))
            f1d.suptitle("\nRadial distribution of scattered photons in detector plane", fontsize=16)
            str_scaling = "binned pixel"
            f1d_ax_left = f1d.add_axes([0.1, 0.1, 0.35, 0.7],title='Radial scaling:' + str_scaling,xlabel="r [" + str_scaling + "]",ylabel="I(r) [photons/" + str_scaling + "]")
            str_scaling = "nyquist pixel"
            f1d_ax_right = f1d.add_axes([0.55, 0.1, 0.35, 0.7],title='Radial scaling:' + str_scaling,xlabel="r [" + str_scaling + "]",ylabel="I(r) [photons/" + str_scaling + "]")
            [plot_args,legend_args] = get_arguments('binned pixel')
            f1d_ax_left.semilogy(*plot_args)
            f1d_ax_left.legend(*legend_args)
            [plot_args,legend_args] = get_arguments('nyquist pixel')
            f1d_ax_right.semilogy(*plot_args)
            f1d_ax_right.legend(*legend_args)
            f1d.show()
            return
        elif scaling == "binned pixel":
            str_scaling = "binned pixel"
            r = numpy.arange(0,len(Ir_sum),1)
        elif scaling == "nyquist pixel":
            str_scaling == "nyquist pixel"
            r = numpy.arange(0,min([self.nyquistpixel_number_x,self.nyquistpixel_number_y])/2,min([self.nyquistpixel_number_x,self.nyquistpixel_number_y])/2/len(Ir_sum))
        elif scaling == "meter":
            str_scaling = "meter"
            r = numpy.arange(0,min([self.pixel_number_x,self.pixel_number_y])/2*self.pixel_size,min([self.pixel_number_x,self.pixel_number_y])/2*self.pixel_size/len(Ir_sum))
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

    def get_nyquist_pixel_size(self):
        return tools.get_nyquist_pixel_size(self.input_object.detector.distance,self.input_object.source.photon.get_wavelength(),self.input_object.sample.get_area())

    def _get_gapsize(self,X_min,X_max,Y_min,Y_max):
        """
        Returns gapsize of pattern in pixels (binned)
        """
        gapsize = 0
        M = self.input_object.detector.mask
        for i in pylab.arange(X_min,X_max+1,1):
            if (M[:,i]==0).all():
                for j in pylab.arange(i,X_max+1,1):
                    if (M[:,j]==1).any():
                        gapsize = j-i
                        break
                break
        for i in pylab.arange(Y_min,Y_max+1,1):
            if (M[i,:]==0).all():
                for j in pylab.arange(i,Y_max+1,1):
                    if (M[j,:]==1).any():
                        gapsize = j-i
                        break
                break
        return gapsize

    def _get_pattern_limits(self):
        """
        Returns spatial limits of pattern in pixels (binned)
        """
        X_min = 0
        Y_min = 0
        X_max = self.amplitudes.shape[1]
        Y_max = self.amplitudes.shape[0]
        M = self.input_object.detector.mask
        for i in pylab.arange(0,M.shape[1],1):
            if (M[:,i]==1).any():
                X_min = i
                break
        for i in M.shape[1]-pylab.arange(1,M.shape[1],1):
            if (M[:,i]==1).any():
                X_max = i
                break
        for i in pylab.arange(0,M.shape[0],1):
            if (M[i,:]==1).any():
                Y_min = i
                break
        for i in M.shape[0]-pylab.arange(1,M.shape[0],1):
            if (M[i,:]==1).any():
                Y_max = i
                break
        return [X_min,X_max,Y_min,Y_max]

    def plot_pattern(self,**kwargs):
        """
        Creates 2-dimensional plot(s) of the distribution of scattered photons.
        Usage: plot_pattern([scaling],[poissonnoise],[logscale],[saturationlevel])
        Keyword arguments:
        - scaling:          'nyquist', 'meter', 'binned pixel' or 'pixel' (default)
        - noise:            'poisson' or 'none' (default)
        - logscale:         False or True (default)
        - saturationlevel:  True or False (default)
        - use_gapmask:      False or True (default)
        """

        scaling = 'pixel'
        scalingargs = ['nyquist','meter','binned pixel','pixel']
        noise = 'none'
        noiseargs = ['poisson','none']
        logscale = True
        logscaleargs = [False,True]
        saturationlevel = False
        saturationlevelargs = [False,True]
        use_gapmask = True
        use_gapmaskargs = [False,True]
        outfile = False
        outfileargs = [True,False]

        I = self.get_intensity_pattern()

        optionkeys = ["scaling","noise","logscale","saturationlevel","use_gapmask","outfile"]
        options = [scaling,noise,logscale,saturationlevel,use_gapmask]
        optionargs = [scalingargs,noiseargs,logscaleargs,saturationlevelargs,use_gapmaskargs,outfileargs]
        keys = kwargs.keys()
        for i in range(0,len(keys)):
            key = keys[i]
            if not key in optionkeys:
                print "ERROR: %s is not a proper key." % key
                return
            keyarg = kwargs[key]
            j = optionkeys.index(key)
            if not keyarg in optionargs[j]:
                print "ERROR: %s is not a proper argument for %s." % (keyarg,key)
                return
            exec "%s = '%s'" % (key,keyarg)
        
        eff_pixel_size_detector = self.input_object.detector.pixel_size * self.input_object.detector.binning
        pixel_size_detector = self.input_object.detector.pixel_size
        pixel_size_nyquist = tools.get_nyquist_pixel_size(self.input_object.detector.distance,self.input_object.source.photon.get_wavelength(),self.input_object.sample.get_area())
        if scaling == "nyquist":
            I *= eff_pixel_size_nyquist**2/pixel_size_detector**2
            u = eff_pixel_size_detector/pixel_size_nyquist
            str_scaling = "Nyquist pixel"
        elif scaling == "meter":
            I /= eff_pixel_size_detector**2
            u = eff_pixel_size_detector
            str_scaling = "m^2"
        elif scaling == "pixel":
            I /= 1.0*self.input_object.detector.binning**2
            u = self.input_object.detector.binning
            str_scaling = scaling
        elif scaling == "binned pixel":
            u = 1.0
            str_scaling = scaling

        I /= u**2

        if noise == "poisson":
            I = abs(pylab.poisson(I))

        if saturationlevel and self.input_object.detector.saturationlevel > 0:
            I *= u**2/(1.0*self.input_object.detector.binning**2)
            I[I>self.input_object.detector.saturationlevel] = self.input_object.detector.saturationlevel
            I /= u**2/(1.0*self.input_object.detector.binning**2)

        [X_min,X_max,Y_min,Y_max] = self._get_pattern_limits()
        xlimit = u*(X_max-X_min)
        ylimit = u*(Y_max-Y_min)
        gapsize = self._get_gapsize(X_min,X_max,Y_min,Y_max)

        I = I[Y_min:Y_max+1,X_min:X_max+1]

        if str_scaling == "binned pixel":
            if self.input_object.detector.binning == 1:
                str_scaling_label = "pixel"
            else:
                str_scaling_label = "%ix%i binned pixel" % (self.input_object.detector.binning,self.input_object.detector.binning)
        else:
            str_scaling_label = str_scaling

        if use_gapmask:
            M = self.input_object.detector.mask.copy()
            if saturationlevel and self.input_object.detector.saturationlevel > 0:
                M[I>=self.input_object.detector.saturationlevel] = 0
            M[M==0] *= pylab.nan

        if logscale:
            I = pylab.log10(I)
            I[I==-pylab.Inf] = I[I!=-pylab.Inf].min()-1.0
            I *= M
            str_Iscaling = r"$\log\left( I \left[ \frac{\mbox{photons}}{\mbox{%s}} \right] \right)$" % str_scaling
        else:
            str_Iscaling = r"$I \left[ \frac{\mbox{photons}}{\mbox{%s}} \right]$" % str_scaling

        Wsizey = 9
        Wsizex = 9
        fsize = 12
        pylab.clf()
        fig = mpy.figure(1,figsize=(Wsizex,Wsizey))
        mpy.rcParams['figure.figsize'] = Wsizex,Wsizey
        fig.suptitle(r"\n - PROPAGATOR -", fontsize=fsize+2)
        alignment = {'horizontalalignment':'center','verticalalignment':'center'}

        fig.text(0.5,(16.75/18.0),r"$E_{\mbox{photon}} = %.0f$ eV ; $\lambda = %.2f$ nm ; $N_{\mbox{photons}} = %.1e$ ; $D_{\mbox{detector}} = %0.3f$ mm" %  (self.input_object.source.photon.get_energy("eV"),self.input_object.source.photon.get_wavelength()/1.0E-09,self.input_object.source.pulse_energy/self.input_object.source.photon.get_energy(),self.input_object.detector.distance/1.0E-03),fontsize=fsize,bbox=dict(fc='0.9',ec="0.9",linewidth=10.0),**alignment) 

        ax = fig.add_axes([3/15.0,5/18.0,10/15.0,10/18.0],title=r'Simulated intensity readout')
        ax.set_xlabel(r"$x$ [" + str_scaling_label + "]",fontsize=fsize)
        ax.set_ylabel(r"$y$ [" + str_scaling_label + "]",fontsize=fsize)

        axcolor = fig.add_axes([3/15.0,3.5/18.0,10/15.0,0.5/18.0])
        for a in [ax,axcolor]:
            for label in a.xaxis.get_ticklabels():
                label.set_fontsize(fsize)
            for label in a.yaxis.get_ticklabels():
                label.set_fontsize(fsize)

        im = ax.matshow(I,extent=[-xlimit/2,xlimit/2,-ylimit/2,ylimit/2],interpolation="nearest",)
        cb = fig.colorbar(im, cax=axcolor,orientation='horizontal')
        cb.set_label(str_Iscaling,fontsize=fsize)

        oversampling_ratio = pixel_size_nyquist/eff_pixel_size_detector
        oversampling_ratio_wo_binning = pixel_size_nyquist/pixel_size_detector
        D = self.input_object.detector.distance
        A =  self.input_object.sample.get_area()
        wavelength = self.input_object.source.photon.get_wavelength()
        [res_horizontally,res_vertically] = self.input_object.get_max_achievable_crystallographic_resolution()
        res_corner = 1/pylab.sqrt(1/res_horizontally**2 + 1/res_vertically**2)
        miss_Ny = gapsize*eff_pixel_size_detector/pixel_size_nyquist
        fig.text(0.5,(1./18.0),r"\textbf{Properties}\\ Linear oversampling ratio: $%.2f$ (binning $%i\times%i$) ; $%.2f$ (no pixel binning)\\" % (oversampling_ratio,self.input_object.detector.binning,self.input_object.detector.binning,oversampling_ratio_wo_binning)+
                 r"Crystallographic resolution (full period): $%.1f$ nm (horizontal) ; $%.1f$ nm (vertical) ; $%.1f$ nm (corner)\\" % (res_horizontally/1.0E-09,res_vertically/1.0E-09,res_corner/1.0E-09)+
                 r"Gap width: $g=%.2f\mbox{ mm}=%.1f$ Nyquist pixels" % (gapsize*eff_pixel_size_detector/1.0E-03,miss_Ny),fontsize=fsize,bbox=dict(fc='0.9',ec="0.9",linewidth=10.0),**alignment)
        #if miss_Ny>2.8:
        #    print "\n!!!\nMissing mode(s) expected (gap width: %.1f Nyquist pixels) \n\nTweaking of one of the parameters recommended:\n- Wavelength w = %.2f nm\n- Sample radius r = %.0f nm\n- Gap size g = %.1f mm\n- Detector distance d = %.0f mm" % (miss_Ny,(rec_wavelength+0.01E-9)*1.0E9,(rec_r-1.0E-9)*1.0E9,(rec_gapsize-0.1E-3)*1.0E3,(rec_d+1.0E-3)*1.0E3)

        if outfile:
            mpy.savefig("intensity_pattern.png",dpi=300)
        else:
            fig.show()
   
    def save_pattern_to_file(self,filename,scaling="binned pixel",*arguments):
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

class Propagation:
    """
    Subclass of the input object.
    Holding propagation parameters.

    """
    
    def __init__(self,**kwargs):
        self._parent = kwargs.get('parent',None)
        self.rs_oversampling = kwargs.get('rs_oversampling',2.0)
        self.N_processes = 4
