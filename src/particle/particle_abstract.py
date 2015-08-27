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

# System packages
import sys, numpy
import copy

# Logging
import logging
logger = logging.getLogger(__name__)
import condor.utils.log
from condor.utils.log import log 
from condor.utils.log import log_and_raise_error,log_warning,log_info,log_debug

# Constants
from scipy import constants

# Condor modules
import condor
from material import Material
from condor.utils.variation import Variation

from condor.utils.config import load_config
import condor.utils.diffraction

class AbstractParticle:
    def __init__(self,
                 rotation_values = None, rotation_formalism = None, rotation_mode = "extrinsic",
                 concentration = 1.,
                 position = None,  position_variation = None, position_spread = None, position_variation_n = None):
        
        self.set_alignment(rotation_values=rotation_values, rotation_formalism=rotation_formalism, rotation_mode=rotation_mode)
        self.set_position_variation(position_variation=position_variation, position_spread=position_spread, position_variation_n=position_variation_n)
        self.position_mean = position if position is not None else [0., 0., 0.]
        self.concentration = concentration
        
    def get_next(self):
        O = {}
        O["_class_instance"]      = self
        O["extrinsic_quaternion"] = self._get_next_extrinsic_rotation().get_as_quaternion()
        O["position"]             = self._get_next_position()
        return O

    def get_current_rotation(self):
        return self._rotations.get_current()

    def set_alignment(self, rotation_values=None, rotation_formalism=None, rotation_mode="extrinsic"):
        # Check input
        if rotation_mode not in ["extrinsic","intrinsic"]:
            log_and_raise_error(logger, "%s is not a valid rotation mode for alignment." % rotation_mode)
            sys.exit(1)
        self._rotation_mode = rotation_mode
        self._rotations = condor.utils.rotation.Rotations(values=rotation_values, formalism=rotation_formalism)

    def set_position_variation(self, position_variation, position_spread, position_variation_n):
        self._position_variation = Variation(position_variation,position_spread,position_variation_n,number_of_dimensions=3,name="particle position")

    
    def _get_next_extrinsic_rotation(self):
        rotation = self._rotations.get_next()
        if self._rotation_mode == "intrinsic":
            rotation = copy.deepcopy(rotation)
            rotation.invert()
        return rotation

    def _get_next_position(self):
        return self._position_variation.get(self.position_mean)
    
    def get_conf(self):
        conf = {}
        conf.update(self._get_conf_rotation())
        conf.update(self._get_conf_position_variation())
        conf["concentration"] = self.concentration
        return conf

    def _get_conf_alignment(self):
        R = self.get_current_rotation()
        A = {
            "rotation_values"    : self._rotations.get_values(),
            "rotation_formalism" : self._rotations.get_formalism(),
            "rotation_mode"      : self._rotation_mode
        }
        return A
    
    def _get_conf_position_variation(self):
        A = {
            "position_variation":        self._position_variation.get_mode(),
            "position_variation_spread": self._position_variation.get_spread(),
            "position_variation_n":      self._position_variation.n
        }
        return A
        

class AbstractContinuousParticle(AbstractParticle):
    def __init__(self,
                 diameter, diameter_variation = None, diameter_spread = None, diameter_variation_n = None,
                 rotation_values = None, rotation_formalism = None, rotation_mode = "extrinsic",
                 concentration = 1.,
                 position = None,  position_variation = None, position_spread = None, position_variation_n = None,
                 material_type = None, massdensity = None, atomic_composition = None):
        
        # Initialise base class
        AbstractParticle.__init__(self,
                                  rotation_values=rotation_values, rotation_formalism=rotation_formalism, rotation_mode=rotation_mode,
                                  concentration=concentration,
                                  position=position, position_variation=position_variation, position_spread=position_spread, position_variation_n=position_variation_n)
        # Diameter
        self.set_diameter_variation(diameter_variation=diameter_variation, diameter_spread=diameter_spread, diameter_variation_n=diameter_variation_n)
        self.diameter_mean = diameter
        # Material
        self.set_material(material_type=material_type, massdensity=massdensity, atomic_composition=atomic_composition)

    def get_conf(self):
        conf = {}
        conf.update(AbstractParticle.get_conf(self))
        conf["diameter"] = self.diameter_mean
        dvar = self._diameter_variation.get_conf()
        conf["diameter_variation"] = dvar["mode"]
        conf["diameter_spread"] = dvar["spread"]
        conf["diameter_variation_n"] = dvar["n"]
        return conf
        
    def get_next(self):
        O = AbstractParticle.get_next(self)
        O["diameter"] = self._get_next_diameter()
        return O

    def set_diameter_variation(self, diameter_variation=None, diameter_spread=None, diameter_variation_n=None):
        self._diameter_variation = Variation(diameter_variation, diameter_spread, diameter_variation_n, name="sample diameter")       

    def _get_next_diameter(self):
        d = self._diameter_variation.get(self.diameter_mean)
        # Non-random diameter
        if self._diameter_variation._mode in [None,"range"]:
            if d <= 0:
                log_and_raise_error(logger,"Sample diameter smaller-equals zero. Change your configuration.")
            else:
                return d
        # Random diameter
        else:
            if d <= 0.:
                log_warning(logger, "Sample diameter smaller-equals zero. Try again.")
                return self._get_next_diameter()
            else:
                return d

    def set_material(self, material_type = None, massdensity = None, atomic_composition = None):
        self.material = Material(material_type=material_type, massdensity=massdensity, atomic_composition=atomic_composition)


