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
from scipy import constants
from python_tools.imgtools import array_to_array      

from variation import Variation
from sample import AbstractSample

class SampleSpheroid(AbstractSample):

    def __init__(self,**kwargs):
        AbstractSample.__init__(self,**kwargs)
        reqk = ["flattening"]
        for k in reqk:
            if k not in kwargs.keys():
                logger.error("Cannot initialize SampleSpheroid instance. %s is a necessary keyword." % k)
                return
        self._flattening_mean = kwargs["flattening"]
        self.set_flattening_variation(kwargs.get("flattening_variation",None),
                                      kwargs.get("flattening_spread",None),
                                      kwargs.get("flattening_variation_n",None))
        self._after_init(**kwargs)

    # Overload
    def _next(self):
        self._next_flattening()
        self._next0()

    def set_flattening_variation(self,flattening_variation=None,flattening_spread=None,flattening_variation_n=None,**kwargs):
        self._flattening_variation = Variation(flattening_variation,flattening_spread,flattening_variation_n,name="spheroid flattening")       

    def _next_flattening(self):
        f = self._flattening_variation.get(self._flattening_mean)
        # Non-random 
        if self._flattening_variation._mode in [None,"range"]:
            if f <= 0:
                logger.error("Spheroid flattening smaller-equals zero. Change your configuration.")
            else:
                self.flattening = f
        # Random 
        else:
            if f <= 0.:
                logger.warning("Spheroid flattening smaller-equals zero. Try again.")
                self._next_flattening()
            else:
                self.flattening = f

    def get_a(self):
        return condortools.to_spheroid_semi_diameter_a(self.diameter,self.flattening)
        
    def get_c(self):
        return condortools.to_spheroid_semi_diameter_c(self.diameter,self.flattening)

    def get_theta(self):
        v_z = numpy.array([1.0,0.0,0.0])
        v_y = numpy.array([0.0,1.0,0.0])
        v_rot = condortools.rotation(v_y,self.euler_angle_0,self.euler_angle_1,self.euler_angle_2)
        theta = numpy.arcsin(numpy.dot(v_rot,v_z))
        return theta

    def get_phi(self):
        v_y = numpy.array([0.0,1.0,0.0])
        v_rot = condortools.rotation(v_y,self.euler_angle_0,self.euler_angle_1,self.euler_angle_2)
        v_rot[0] = 0.0
        v_rot = v_rot / numpy.sqrt(v_rot[0]**2+v_rot[1]**2+v_rot[2]**2)       
        phi = numpy.arccos(numpy.dot(v_rot,v_y))
        return phi

    def propagate_single(self,detector0=None,source0=None):
        # scattering amplitude from homogeneous sphere
        if source0 == None:
            source = self._parent.source
        else:
            source = source0
        if detector0 == None:
            detector = self._parent.detector
        else:
            detector = detector0

        dn = self._get_dn()
        F0 = self._get_F0(source,detector)
        # Rotation is taken into account directly in the diffraction formula (much faster)
        q = detector.generate_qmap(euler_angle_0=0.,euler_angle_1=0.,euler_angle_2=0.)
        qx = q[:,:,2]
        qy = q[:,:,1]

        R = self.diameter/2.
        V = 4/3.*numpy.pi*R**3
        K = (F0*V*dn.real)**2
        F = condortools.F_spheroid_diffraction(K,qx,qy,self.get_a(),self.get_c(),self.get_theta(),self.get_phi())

        return {"amplitudes":F,"sample_diameter":self.diameter,"euler_angle_0":self.euler_angle_0,"euler_angle_1":self.euler_angle_1,"euler_angle_2":self.euler_angle_2,
                "sample_spheroid_phi":self.get_phi(),"sample_spheroid_theta":self.get_theta(),"sample_spheroid_flattening":self.flattening,"sample_spheroid_diameter_a":self.get_a()*2,"sample_spheroid_diameter_c":self.get_c()*2}

    def get_area(self):
        """
        Calculates area of projected spheroid
        """
        logger.warning("Calculates area of WRONGLY projected spheroid, fix when there is time.")
        return ((4/3.*numpy.pi*self.get_a()**2*self.get_c())**(2/3.))
