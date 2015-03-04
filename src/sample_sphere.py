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

from sample import AbstractSample

class SampleSphere(AbstractSample):
    """
    A class of the input-object.
    Sample is a homogeneous sphere defined by its diameter and a material object.

    """

    def __init__(self,**kwargs):
        AbstractSample.__init__(self,**kwargs)
        self._after_init(**kwargs)

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
        q = detector.generate_absqmap()
        
        R = self.diameter/2.
        V = 4/3.*numpy.pi*R**3
        K = (F0*V*dn.real)**2
        #K = source.get_intensity()*(self.material.get_electron_density()*detector.get_pixel_size("binned")/detector.distance*constants.value("classical electron radius")*V)**2
        F = condortools.F_sphere_diffraction(K,q,R)

        return {"amplitudes":F,"sample_diameter":self.diameter}

    def get_area(self):
        """ Calculates area of projected sphere """
        return numpy.pi*(self.diameter/2.)**2
