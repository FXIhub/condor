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
        self._theta_mean      = kwargs.get("theta",0.)
        self._phi_mean        = kwargs.get("phi",0.)
        self.set_theta_variation(kwargs.get("theta_variation",None),
                                 kwargs.get("theta_spread",None),
                                 kwargs.get("theta_variation_n",None))
        self.set_phi_variation(kwargs.get("phi_variation",None),
                                 kwargs.get("phi_spread",None),
                                 kwargs.get("phi_variation_n",None))
        self._after_init(**kwargs)

    # Overload
    def _next(self):
        self._next_flattening()
        self._next_theta()
        self._next_phi()
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

    def _set_geometry(self,diameter,flattening):
        self._a = condortools.to_spheroid_semi_diameter_a(diameter,flattening)
        self._c = condortools.to_spheroid_semi_diameter_c(diameter,flattening)

    def set_theta_variation(self,theta_variation=None,theta_spread=None,theta_variation_n=None,**kwargs):
        self._theta_variation = Variation(theta_variation,theta_spread,theta_variation_n,name="spheroid theta")       

    def _next_theta(self):
        self.theta = self._theta_variation.get(self._theta_mean)

    def set_phi_variation(self,phi_variation=None,phi_spread=None,phi_variation_n=None,**kwargs):
        self._phi_variation = Variation(phi_variation,phi_spread,phi_variation_n,name="spheroid phi")       

    def _next_phi(self):
        self.phi = self._phi_variation.get(self._phi_mean)

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

        V = 4/3.*numpy.pi*self.a**2*self.c
        dn = self._get_dn()
        F0 = self._get_F0(source,detector)
        K = (F0*V*dn.real)**2
        
        q = detector.generate_qmap(euler_angle_0=self.euler_angle_0,euler_angle_1=self.euler_angle_1,euler_angle_2=self.euler_angle_2)
        qx = q[:,:,2]
        qy = q[:,:,1]
        F = condortools.F_spheroid_diffraction(K,qx,qy,self.a,self.c,self.theta,self.phi)

        return {"amplitudes":F,"euler_angle_0":self.euler_angle_0,"euler_angle_1":self.euler_angle_1,"euler_angle_2":euler_angle_2}

    def get_area(self):
        """
        Calculates area of projected spheroid
        """
        logger.warning("Calculates area of WRONGLY projected spheroid, fix when there is time.")
        return ((4/3.*numpy.pi*self.a**2*self.c)**(2/3.))
