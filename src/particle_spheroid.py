# -----------------------------------------------------------------------------------------------------
# CONDOR
# Simulator for diffractive single-particle imaging experiments with X-ray lasers
# http://xfel.icm.uu.se/condor/
# -----------------------------------------------------------------------------------------------------
# Copyright 2014 Max Hantke, Filipe R.N.C. Maia, Tomas Ekeberg
# Condor is distributed under the terms of the GNU General Public License
# -----------------------------------------------------------------------------------------------------
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but without any warranty; without even the implied warranty of
# merchantability or fitness for a pariticular purpose. See the
# GNU General Public License for more details.
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
# -----------------------------------------------------------------------------------------------------
# General note:
# All variables are in SI units by default. Exceptions explicit by variable name.
# -----------------------------------------------------------------------------------------------------

import numpy

import logging
logger = logging.getLogger("Condor")
import utils.log
from utils.log import log 

from particle_abstract import AbstractContinuousParticleModel

from utils.variation import Variation


class ParticleModelSpheroid(AbstractContinuousParticleModel):
    def __init__(self,**kwargs):
        # Check for valid set of keyword arguments
        self.req_keys += ["flattening"]
        self.opt_keys += ["flattening_variation","flattening_spread","flattening_variation_n"]
        # Start initialisation
        AbstractContinuousParticleModel.__init__(self,**kwargs)
        self.flattening_mean = kwargs["flattening"]
        self.set_flattening_variation(flattening_variation=kwargs.get("flattening_variation",None),flattening_spread=kwargs.get("flattening_spread",None),flattening_variation_n=kwargs.get("flattening_variation_n",None))

    def get_next(self):
        O = AbstractContinuousParticleModel.get_next(self)
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
