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

import sys,numpy
import logging
logger = logging.getLogger("Condor")
if "utils" not in sys.path: sys.path.append("utils")
import condortools
from variation import Variation
from material import Material

class AbstractSample:
    def __init__(self,**kwargs):
        self._parent = kwargs.get('parent',None)
        # Diameter
        reqk = ["diameter"]
        for k in reqk:
            if k not in kwargs.keys():
                logger.error("Cannot initialize AbstractSample instance. %s is a necessary keyword." % k)
                return
        self._diameter_mean = kwargs['diameter']
        # Material
        materialargs = {}
        if 'massdensity' in kwargs:
            materialargs['massdensity'] = kwargs['massdensity']
            for key in kwargs.keys():
                if key[0] == 'c' and len(key) <= 3: materialargs[key] = kwargs[key]
            self.material = Material(self,**materialargs)
        elif "material_type" in kwargs:
            materialargs['material_type'] = kwargs['material_type']
            self.material = Material(self,**materialargs)
        else:
            self.material = None
        # Number of images
        if "number_of_images" in kwargs:
            self.number_of_images = kwargs["number_of_images"]
        else:
            # Maintaining depreciated keyword
            self.number_of_images = kwargs.get("number_of_orientations",1)
        self._after_init_called = False

    # This has to be called by every inheriting class after init!
    def _after_init(self,**kwargs):
        self.set_alignment(**kwargs)
        self.set_diameter_variation(**kwargs)
        self._i = -1
        self._next()
        self._after_init_called = True

    def propagate(self,detector0=None,source0=None,output=["amplitudes"]):
        if not self._after_init_called:
            logger.error("Sample class was not successfully initialised. Cannot perform propagation.")
            exit(1)

        if source0 is None:
            source = self._parent.source
        else:
            source = source0
        if detector0 is None:
            detector = self._parent.detector
        else:
            detector = detector0

        O_all = {}
        while self._i < self.number_of_images:
            logger.info("Calculation diffraction pattern (%i/%i). (PROGSTAT)" % (self._i+1,self.number_of_images))
            O = self.propagate_single()
            # Collect values of sample variables
            for k,d in O.items():
                if k not in O_all:
                    O_all[k] = []
                O_all[k].append(d)

            # Collect values of detector variables
            for k in ["cx","cy","cxXxX","cyXxX"]:
                if k not in O_all:
                    O_all[k] = []
            O_all["cx"].append(detector.cx)
            O_all["cy"].append(detector.cy)
            O_all["cxXxX"].append(detector.get_cx("binned"))
            O_all["cyXxX"].append(detector.get_cy("binned"))
            # Collect values of source variables
            for k in ["beam_intensity","beam_intensity_ph_per_um2"]:
                if k not in O_all:
                    O_all[k] = []
            O_all["beam_intensity_ph_per_um2"].append(source.get_intensity("ph/m2"))
            O_all["beam_intensity"].append(source.get_intensity("J/m2"))
            # Iteration step
            self._next()
            detector._next()
            source._next()
        for k,i in O_all.items():
            O_all[k] = numpy.array(O_all[k])
        return O_all

    # Overload this function!
    #def propagate_single(self,detector0=None,source0=None):
    #    pass

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

    def set_alignment(self,alignment=None,**kwargs):
        if alignment not in [None,"random","euler_angles","first_axis"]:
            logger.error("Invalid argument for sample alignment specified.")
            return
        self._euler_angle_0 = kwargs.get("euler_angle_0",None)
        self._euler_angle_1 = kwargs.get("euler_angle_1",None)
        self._euler_angle_2 = kwargs.get("euler_angle_2",None)
        if alignment is None and (self._euler_angle_0 is not None and self._euler_angle_1 is not None and self._euler_angle_2 is not None):
            print "bla"
            self._alignment = "euler_angles"
        else:
            self._alignment = alignment
        
    # Overload this function if needed
    def _next(self):
        self._next0()

    def _next0(self):
        self._next_orientation()
        self._next_diameter()
        self._i += 1

    def _next_orientation(self):
        if self._alignment == "first_axis":
            # Sanity check
            if self._euler_angle_0 is not None or self._euler_angle_1 is not None or self._euler_angle_2 is not None:
                logger.error("Conflict of arguments: Specified first_axis alignment and also specified set of euler angles. This does not make sense.")
                return
            self.euler_angle_0 = 0.
            self.euler_angle_1 = 0.
            self.euler_angle_2 = 0.
        elif self._alignment == "random":
            # Sanity check
            if self._euler_angle_0 is not None or self._euler_angle_1 is not None or self._euler_angle_2 is not None:
                logger.error("Conflict of arguments: Specified random alignment and also specified set of euler angles. This does not make sense.")
                return
            (self.euler_angle_0,self.euler_angle_1,self.euler_angle_2) = condortools.random_euler_angles()
        elif self._alignment == "euler_angles":
            # Many orientations (lists of euler angles)
            if isinstance(self._euler_angle_0,list):
                self.euler_angle_0 = self._euler_angle_0[self._i]
                self.euler_angle_1 = self._euler_angle_1[self._i]
                self.euler_angle_2 = self._euler_angle_2[self._i]
            # One orientation (euler angles are scalars)
            else:
                self.euler_angle_0 = self._euler_angle_0
                self.euler_angle_1 = self._euler_angle_1
                self.euler_angle_2 = self._euler_angle_2

    def set_diameter_variation(self,diameter_variation=None,diameter_spread=None,diameter_variation_n=None,**kwargs):
        self._diameter_variation = Variation(diameter_variation,diameter_spread,diameter_variation_n,name="sample diameter")       

    def _next_diameter(self):
        d = self._diameter_variation.get(self._diameter_mean)
        # Non-random diameter
        if self._diameter_variation._mode in [None,"range"]:
            if d <= 0:
                logger.error("Sample diameter smaller-equals zero. Change your configuration.")
            else:
                self.diameter = d
        # Random diameter
        else:
            if d <= 0.:
                logger.warning("Sample diameter smaller-equals zero. Try again.")
                self._next_diameter()
            else:
                self.diameter = d




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
