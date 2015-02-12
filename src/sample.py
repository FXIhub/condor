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

import sys,numpy,time,multiprocessing
import logging
logger = logging.getLogger("Condor")
if "utils" not in sys.path: sys.path.append("utils")
import config,imgutils,condortools
import utils.nfft
import utils.icosahedron

# Pythontools
from python_tools import gentools,cxitools,imgtools


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
                    exec "self." + key + " = kwargs[key]"
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
            self.cAu = config.DICT_atomic_composition[self.material_type][6]
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

    def _get_dn(self):
        # refractive index from material
        if self.material == None:
            dn = complex(1.,0.)
        else:
            dn = self.material.get_dn()
        return dn

    def _get_F0(self,source0=None,detector0=None):
        if source0 == None:
            source = self._parent.source
        else:
            source = source0
        if detector0 == None:
            detector = self._parent.detector
        else:
            detector = detector0
        # F0 = sqrt(I_0 Omega_p) 2pi/wavelength^2
        wavelength = source.photon.get_wavelength()
        I_0 = source.get_intensity("ph/m2")
        Omega_p = detector.get_pixel_solid_angle("binned")
        F0 = numpy.sqrt(I_0*Omega_p)*2*numpy.pi/wavelength**2 
        return F0

    def set_random_orientation(self):
        self.set_alignment("random")

    def set_alignment(self,alignment="first_axis",euler_angle_0=None,euler_angle_1=None,euler_angle_2=None):
        if alignment not in ["random","euler_angles","first_axis"]:
            logger.error("Invalid argument for sample alignment specified.")
            return
        self._alignment = kwargs["alignment"]
        self._euler_angle_0 = euler_angle_0
        self._euler_angle_1 = euler_angle_1
        self._euler_angle_2 = euler_angle_2

    def set_diameter_variation(self,variation="none",spread=0):
        if variation not in ["none","uniform","normal"]:
            logger.error("Invalid argument for diameter variation specified.")
            return
        self._diameter_variation = variation
        self._diameter_spread = 0.

    def _get_diameters(self,N=None):
        if N is None:
            _N = 1
        else:
            _N = N
        if self._diameter_variation == "none":
            d = list(numpy.ones(_N) * self.radius * 2)
        elif self._diameter_variation == "uniform":
            d = list(numpy.random.uniform(self.radius*2-self._diameter_spread/2.,self.radius*2+self._diameter_spread/2.,_N))
        elif self._diameter_variation == "normal":
            d = list(numpy.random.normal(self.radius*2,self._diameter_spread,_N))
        return d

    def _get_euler_angles(self,N=None):        
        if self._alignment == "first_axis":
            if N is None:
                _N = 1
            else:
                _N = N
            e0 = list(numpy.zeros(_N))
            e1 = list(numpy.zeros(_N))
            e2 = list(numpy.zeros(_N))
        elif self._alignment == "random":
            # Sanity check
            if self._euler_angle_0 != 0. or self._euler_angle_1 != 0. or self._euler_angle_2 != 0.:
                logger.error("Conflict of arguments: Specified random alignment and also specified set of euler angles. This does not make sense.")
                return
            if N is None:
                _N = 1
            else:
                _N = N
            e0 = numpy.zeros(_N)
            e1 = numpy.zeros(_N)
            e2 = numpy.zeros(_N)
            for i in range(_N):
                (e0[i],e1[i],e2[i]) = condortools.random_euler_angles()
            e0 = list(e0)
            e1 = list(e1)
            e2 = list(e2)
        elif self._alignment == "euler_angles":
            # Many orientations (lists of euler angles)
            if isinstance(euler_angle_0,list):
                # Sanity check
                if N is not None:
                    if len(euler_angle_0) != N:
                        logger.error("Conflict of arguments: N = %i and len(euler_angle_0) = %i." % (N,len(euler_angle_0)))
                        return
                _N = len(self._euler_angle_0)
                e0 = self._euler_angle_0
                e1 = self._euler_angle_1
                e2 = self._euler_angle_2
            # One orientation (euler angles are scalars)
            else:
                _N = 1
                e0 = [self._euler_angle_0]
                e1 = [self._euler_angle_1]
                e2 = [self._euler_angle_2]
        return (e0,e1,e2)

        


class SampleSphere(Sample):
    """
    A class of the input-object.
    Sample is a homogeneous sphere defined by a radius and a material object.

    """

    def __init__(self,**kwargs):
        Sample.__init__(self,**kwargs)
        reqk = ["diameter"]
        for k in reqk:
            if k not in kwargs.keys():
                logger.error("Cannot initialize SampleSphere instance. %s is a necessary keyword." % k)
                return
        self.radius = kwargs['diameter']/2.
        self._parent = kwargs.get('parent',None)

        material_kwargs = kwargs.copy()
        non_material_arguments = ['parent','diameter','diameter_variation','diameter_spread']
        for non_material_argument in non_material_arguments:
            if non_material_argument in material_kwargs:
                material_kwargs.pop(non_material_argument)
        material_kwargs['parent'] = self
        self.material = Material(**material_kwargs)

    def propagate(self,detector0=None,source0=None,number_of_images=None):
        # scattering amplitude from homogeneous sphere
        if source0 == None:
            source = self._parent.source
        else:
            source = source0
        if detector0 == None:
            detector = self._parent.detector
        else:
            detector = detector0

        F = []
        dn = self._get_dn()
        F0 = self._get_F0(source,detector)
        q = detector.generate_absqmap()

        d = self.get_diameters(number_of_images)

        for i in range(len(d)):
            R = d[i]/2.
            V = 4/3.*numpy.pi*R**3
            K = (F0*V*dn.real)**2
            F.append(condortools.F_sphere_diffraction(K,q,R))

        return {"amplitudes":F,"diameters":d}

    def get_area(self):
        """ Calculates area of projected sphere """
        return numpy.pi*self.radius**2

class SampleSpheroid(Sample):

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

    def propagate(self,detector0=None,source0=None,number_of_images=None):
        # scattering amplitude from homogeneous sphere
        if source0 == None:
            source = self._parent.source
        else:
            source = source0
        if detector0 == None:
            detector = self._parent.detector
        else:
            detector = detector0

        V = 4/3.*numpy.pi*self.a**2*self.c
        dn = self._get_dn()
        F0 = self._get_F0(source,detector)
        K = (F0*V*dn.real)**2
        q = detector.generate_qmap(euler_angle_0=0.,euler_angle_1=0.,euler_angle_2=0.)
        qx = q[:,:,2]
        qy = q[:,:,1]
        F = [condortools.F_spheroid_diffraction(K,qx,qy,self.a,self.c,self.theta,self.phi)]

        return {"amplitudes":F}

    def get_area(self):
        """
        Calculates area of projected spheroid
        """
        logger.warning("Calculates area of WRONGLY projected spheroid, fix when there is time.")
        return (4/3.*numpy.pi*self.a**2*self.c)**(2/3.)


class SampleMap(Sample):

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
          - diameter_c: diameter along the singular axis (rotation axis of ellipsoid)
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
        self.map3d_fine = None
        self._map3d = None
        self._dX = None
        self._map3d_fine = None
        self._dX_fine = None
        self.radius = kwargs.get('radius',None)

        if "dx_fine" in kwargs:
            self.dX_fine = kwargs["dx_fine"]
        elif "oversampling_fine" in kwargs:
            #self.dX_fine = self._parent.detector.get_real_space_resolution_element()/float(kwargs["oversampling_fine"])
            self.dX_fine = self._parent.detector.get_real_space_resolution_element()/float(kwargs["oversampling_fine"])/numpy.sqrt(2)

        # Map
        if "geometry" in kwargs:
            if kwargs["geometry"] == "icosahedron":
                if "diameter" not in kwargs:
                    logger.error("Cannot initialize SampleMap instance. diameter is a necessary keyword for geometry=icosahedron.") 
                self.put_icosahedron(kwargs["diameter"]/2.,**kwargs)
                self.radius = kwargs["diameter"]/2.
            elif kwargs["geometry"] == "spheroid":
                if "diameter_a" not in kwargs or "diameter_c" not in kwargs:
                    logger.error("Cannot initialize SampleMap instance. a_diameter and c_diameter are necessary keywords for geometry=spheroid.")
                self.put_spheroid(kwargs["diameter_a"]/2.,kwargs["diameter_c"]/2.,**kwargs)
                self.radius = (2*kwargs["diameter_a"]+kwargs["diameter_c"])/3./2.
            elif kwargs["geometry"] == "sphere":
                if "diameter" not in kwargs:
                    logger.error("Cannot initialize SampleMap instance. diameter is a necessary keyword for geometry=sphere.")
                self.put_sphere(kwargs["diameter"]/2.,**kwargs)
                self.radius = kwargs["diameter"]/2.
            if kwargs["geometry"] == "cube":
                if "edge_length" not in kwargs:
                    logger.error("Cannot initialize SampleMap instance. edge_length is a necessary keyword for geometry=cube.") 
                self.put_cube(kwargs["edge_length"],**kwargs)
                self.radius = (4/3/numpy.pi)**(1/3.)*kwargs["edge_length"] # volume equivalent radius
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
                    self.dX_fine = kwargs["diameter"]/float(s[0])
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

        self.set_alignment(**kwargs)

    def propagate(self,detector0=None,source0=None,number_of_images=None):
        # scattering amplitude from dn-map: F = F0 DFT{dn} dV
        if source0 == None:
            source = self._parent.source
        else:
            source = source0
        if detector0 == None:
            detector = self._parent.detector
        else:
            detector = detector0

        map3d = None
        #dX = detector.get_real_space_resolution_element()
        dX = detector.get_real_space_resolution_element() / numpy.sqrt(2)
        self.dX = dX

        if self.dX_fine > dX:
            logger.error("Finer real space sampling required for chosen geometry.")
            return

        # has map3d_fine the required real space grid?
        if map3d == None and abs(self.dX_fine/dX-1) < 0.001:
            # ok, we'll take the fine map
            map3d = self.map3d_fine
            logger.debug("Using the fine map for propagtion.")
            self._map3d = self.map3d_fine

        # do we have an interpolated map?
        if map3d == None and self._dX != None:
            # does it have the right spacing?
            if abs(self._dX/dX-1) < 0.001:
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
        if map3d == None and self.dX_fine < dX:
            from scipy import ndimage
            f = self.dX_fine/dX
            N_mapfine = self.map3d_fine.shape[0]
            L_mapfine = (N_mapfine-1)*self.dX_fine
            N_map = int(numpy.floor((N_mapfine-1)*f))+1
            L_map = (N_map-1)*dX
            gt = numpy.float64(numpy.indices((N_map,N_map,N_map)))/float(N_map-1)*(N_mapfine-1)*L_map/L_mapfine
            map3d = ndimage.map_coordinates(self.map3d_fine, gt, order=3)
            # Cache interpolated data 
            self._map3d = map3d
            self._dX = N_mapfine/(1.*N_map)*self.dX_fine
            # Cace fine data for later decision whether or not the interpolated map can be used again
            self._map3d_fine = self.map3d_fine
            self._dX_fine = self.dX_fine
            logger.debug("Using a newly interpolated map for propagtion.")

        dn_map3d = numpy.array(map3d,dtype="complex128") * self._get_dn()
        self.dn_map3d = dn_map3d

        (e0,e1,e2) = self._get_euler_angles(number_of_images)

        F = []
        for i in range(len(e0)):
            logger.info("Calculation diffraction pattern (%i/%i). (PROGSTAT)" % (i+1,len(e0)))

            # scattering vector grid
            q_scaled = detector.generate_qmap(nfft_scaled=True,euler_angle_0=e0[i],euler_angle_1=e1[i],euler_angle_2=e2[i])
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
            F.append(self._get_F0(source,detector) * fourierpattern * dX**3)
    
        return {"amplitudes": F, "euler_angle_0": e0[i], "euler_angle_1": e1[i], "euler_angle_2": e2[i], "F0": self._get_F0(source, detector) , "dX3": dX**3, "grid": q_reshaped, 'qmap3d': qmap3d}
        
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
        #print map_add.shape,origin
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
 
    def put_spheroid(self,a,b,**kwargs):
        e0 = kwargs.get("geometry_euler_angle_0")
        e1 = kwargs.get("geometry_euler_angle_1")
        e2 = kwargs.get("geometry_euler_angle_2")
        # maximum radius
        Rmax = max([a,b])
        # maximum radius in pixel
        nRmax = Rmax/self.dX_fine
        # dimensions in pixel
        nA = a/self.dX_fine
        nB = b/self.dX_fine
        # leaving a bit of free space around spheroid
        N = int(round((nRmax*1.2)*2))
        spheromap = make_spheroid_map(N,nA,nB,e0,e1,e2)
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
        if self.radius != None: return numpy.pi*self.radius**2

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
    e_c = condortools.rotation(numpy.array([0.0,0.0,1.0]),euler0,euler1,euler2)
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


#def calculatePattern_SampleSpheres(input_obj):
#    wavelength = input_obj.source.photon.get_wavelength()
#    I_0 = input_obj.source.get_intensity("ph/m2")
#    Omega_p = input_obj.detector.get_pixel_solid_angle("binned")
#    if self.material == None:
#        dn_real = 1.
#    else:
#        dn_real = self.material.get_dn().real

#    radii = self.get_radii()
#    dn_real = (1-self.material.get_n()).real
#    absq = input_obj.detector.generate_absqmap()
#    q = input_obj.detector.generate_qmap()
#    F = numpy.zeros(shape=absq.shape,dtype='complex')
#    for R in radii:
#        Fr = numpy.zeros_like(F)
#        Fr[absq!=0.0] = (numpy.sqrt(I_0*Omega_p)*2*numpy.pi/wavelength**2*4/3.0*numpy.pi*R**3*3*(numpy.sin(absq[absq!=0.0]*R)-absq[absq!=0.0]*R*numpy.cos(absq[absq!=0.0]*R))/(absq[absq!=0.0]*R)**3*dn_real)
#        try: Fr[absq==0] = numpy.sqrt(I_0*Omega_p)*2*numpy.pi/wavelength**2*4/3.0*numpy.pi*R**3*dn_real
#        except: pass

#        indices = self.r==R

#        for i in range(sum(indices)):
#            looger.debug("%i" % i)
#            d = [numpy.array(self.z)[indices][i],
#                 numpy.array(self.y)[indices][i],
#                 numpy.array(self.x)[indices][i]]
#            F[:,:] += (Fr*input_obj.get_phase_ramp(q,d))[:,:]
#    return F
