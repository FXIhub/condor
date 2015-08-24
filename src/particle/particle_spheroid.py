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
logger = logging.getLogger(__name__)

import condor
import condor.utils.log
from condor.utils.log import log_and_raise_error,log_warning,log_info,log_debug

from particle_abstract import AbstractContinuousParticle

from condor.utils.variation import Variation


class ParticleSpheroid(AbstractContinuousParticle):
    def __init__(self,
                 diameter,
                 diameter_variation = None, diameter_spread = None, diameter_variation_n = None,
                 flattening = 1., flattening_variation = None, flattening_spread = None, flattening_variation_n = None,
                 rotation_values = None, rotation_formalism = None, rotation_mode = "extrinsic",
                 concentration = 1.,
                 position = None, position_variation = None, position_spread = None, position_variation_n = None,
                 material_type = None, massdensity = None, atomic_composition = None):

        # Initialise base class
        AbstractContinuousParticle.__init__(self,
                                            diameter=diameter, diameter_variation=diameter_variation, diameter_spread=diameter_spread, diameter_variation_n=diameter_variation_n,
                                            rotation_values=rotation_values, rotation_formalism=rotation_formalism, rotation_mode=rotation_mode,
                                            concentration=concentration,
                                            position=position, position_variation=position_variation, position_spread=position_spread, position_variation_n=position_variation_n,
                                            material_type=material_type, massdensity=massdensity, atomic_composition=atomic_composition)
        self.flattening_mean = flattening
        self.set_flattening_variation(flattening_variation=flattening_variation, flattening_spread=flattening_spread, flattening_variation_n=flattening_variation_n)

    def get_conf(self):
        conf = {}
        conf.update(AbstractContinuousParticle.get_conf(self))
        conf["flattening"] = self.flattening_mean
        fvar = self._flattening_variation.get_conf()
        conf["flattening_variation"] = fvar["mode"]
        conf["flattening_spread"] = fvar["spread"]
        conf["flattening_variation_n"] = fvar["n"]
        return conf
        
    def get_next(self):
        O = AbstractContinuousParticle.get_next(self)
        O["flattening"] = self._get_next_flattening()
        return O
        
    def set_flattening_variation(self, flattening_variation=None, flattening_spread=None, flattening_variation_n=None):
        self._flattening_variation = Variation(flattening_variation, flattening_spread, flattening_variation_n, name="spheroid flattening")       

    def _get_next_flattening(self):
        f = self._flattening_variation.get(self.flattening_mean)
        # Non-random 
        if self._flattening_variation._mode in [None, "range"]:
            if f <= 0:
                log_and_raise_error(logger, "Spheroid flattening smaller-equals zero. Change your configuration.")
            else:
                return f
        # Random 
        else:
            if f <= 0.:
                log_warning(logger, "Spheroid flattening smaller-equals zero. Try again.")
                return self._get_next_flattening()
            else:
                return f
