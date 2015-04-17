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

import numpy

import logging
logger = logging.getLogger("Condor")
import utils.log
from utils.log import log 

from particle_abstract import AbstractContinuousParticleSpecies

from utils.variation import Variation


class ParticleSpeciesSpheroid(AbstractContinuousParticleSpecies):
    def __init__(self,**kwargs):
        # Check for valid set of keyword arguments
        self.req_keys = ["flattening"]
        self.opt_keys = ["flattening_variation","flattening_spread","flattening_variation_n"]
        # Start initialisation
        AbstractContinuousParticleSpecies.__init__(self,**kwargs)
        self.flattening_mean = kwargs["flattening"]
        self.set_flattening_variation(flattening_variation=kwargs.get("flattening_variation",None),flattening_spread=kwargs.get("flattening_spread",None),flattening_variation_n=kwargs.get("flattening_variation_n",None))

    def get_next(self):
        O = AbstractContinuousParticleSpecies.get_next(self)
        O["flattening"] = self._get_next_flattening()
        return O
        
    def set_flattening_variation(self,flattening_variation=None,flattening_spread=None,flattening_variation_n=None,**kwargs):
        self._flattening_variation = Variation(flattening_variation,flattening_spread,flattening_variation_n,name="spheroid flattening")       

    def _get_next_flattening(self):
        f = self._flattening_variation.get(self.flattening_mean)
        # Non-random 
        if self._flattening_variation._mode in [None,"range"]:
            if f <= 0:
                log(logger.error,"Spheroid flattening smaller-equals zero. Change your configuration.")
            else:
                return f
        # Random 
        else:
            if f <= 0.:
                log(logger.warning,"Spheroid flattening smaller-equals zero. Try again.")
                return self._get_next_flattening()
            else:
                return f
