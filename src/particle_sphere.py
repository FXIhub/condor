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

class ParticleSpeciesSphere(AbstractContinuousParticleSpecies):
    def __init__(self,**kwargs):
        self.req_keys = []
        self.opt_keys = []
        AbstractContinuousParticleSpecies.__init__(self,**kwargs)
