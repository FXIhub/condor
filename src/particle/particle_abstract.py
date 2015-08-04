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

# Logging
import condor.utils.log
from condor.utils.log import log 

# Constants
from scipy import constants

# Condor modules
import condor
from material import Material
from condor.utils.variation import Variation

import condor.utils.config
import condor.utils.diffraction

def load_particles(conf):
    C = condor.utils.config.load_config(conf)
    names = [k for k in C.keys() if k.startswith("particle")]
    particles = {}
    for n in names:
        particles[n] = load_particle(conf, n)
    return particles

def load_particle(conf, name=None):
    C = condor.utils.config.load_config(conf)
    names = [k for k in C.keys() if k.startswith("particle")]
    default = condor.utils.config.get_default_conf()
    if len(names) > 1 and name is None:
        log(condor.CONDOR_logger.error, "There is more than one particle defined in configuration and no \'name\' is given to decide which one to load.")
        return
    if name is None:
        k = C.keys()[0]
    else:
        k = name
    if k.startswith("particle_sphere"):
        CP = condor.utils.config.load_config({"particle": C[k]}, {"particle": default["particle_sphere"]})
        particle = condor.ParticleSphere(**CP["particle"])
    elif k.startswith("particle_spheroid"):
        CP = condor.utils.config.load_config({"particle": C[k]}, {"particle": default["particle_spheroid"]})
        particle = condor.ParticleSpheroid(**CP["particle"])
    elif k.startswith("particle_map"):
        CP = condor.utils.config.load_config({"particle": C[k]}, {"particle": default["particle_map"]})
        particle = condor.ParticleMap(**CP["particle"])
    elif k.startswith("particle_molecule"):
        CP = condor.utils.config.load_config({"particle": C[k]}, {"particle": default["particle_molecule"]})
        particle = condor.ParticleMolecule(**CP["particle"])
    else:
        log(condor.CONDOR_logger.error,"Particle model for %s is not implemented." % k)
        sys.exit(1)
    return particle

class AbstractParticle:
    def __init__(self,
                 alignment = None, euler_angle_0 = None, euler_angle_1 = None, euler_angle_2 = None,
                 concentration = 1.,
                 position = None,  position_variation = None, position_spread = None, position_variation_n = None):
        
        self.set_alignment(alignment=alignment, euler_angle_0=euler_angle_0, euler_angle_1=euler_angle_1, euler_angle_2=euler_angle_2)
        self.set_position_variation(position_variation=position_variation, position_spread=position_spread, position_variation_n=position_variation_n)
        self.position_mean = position if position is not None else [0., 0., 0.]
        self.concentration = concentration

    def get_conf(self):
        conf = {}
        conf.update(self.get_alignment())
        conf.update(self.get_position_variation())
        conf["concentration"] = self.concentration
        return conf
        
    def get_next(self):
        O = {}
        O["_class_instance"] = self 
        euler_angle_0,euler_angle_1,euler_angle_2 = self._get_next_orientation()
        O["euler_angle_0"] = euler_angle_0
        O["euler_angle_1"] = euler_angle_1
        O["euler_angle_2"] = euler_angle_2
        O["position"] = self._get_next_position()
        return O
    
    def set_alignment(self, alignment=None, euler_angle_0=None, euler_angle_1=None, euler_angle_2=None):
        if alignment not in [None,"random","euler_angles","first_axis","random_euler_angle_0"]:
            log(condor.CONDOR_logger.error,"Invalid argument for sample alignment specified.")
            return
        self._euler_angle_0 = euler_angle_0
        self._euler_angle_1 = euler_angle_1
        self._euler_angle_2 = euler_angle_2
        if alignment is None and (self._euler_angle_0 is not None and self._euler_angle_1 is not None and self._euler_angle_2 is not None):
            self._alignment = "euler_angles"
        else:
            self._alignment = alignment

    def get_alignment(self):
        A = {
            "alignment":     self._alignment,
            "euler_angle_0": self._euler_angle_0,
            "euler_angle_1": self._euler_angle_1,
            "euler_angle_2": self._euler_angle_2,
        }
        return A

    def _get_next_orientation(self):
        if self._alignment == "first_axis" or self._alignment is None:
            # Sanity check
            if self._euler_angle_0 is not None or self._euler_angle_1 is not None or self._euler_angle_2 is not None:
                log(condor.CONDOR_logger.error,"Conflict of arguments: Specified first_axis alignment and also specified set of euler angles. This does not make sense.")
                exit(1)
            euler_angle_0 = 0.
            euler_angle_1 = 0.
            euler_angle_2 = 0.
        elif self._alignment == "random":
            # Sanity check
            if self._euler_angle_0 is not None or self._euler_angle_1 is not None or self._euler_angle_2 is not None:
                log(condor.CONDOR_logger.error,"Conflict of arguments: Specified random alignment and also specified set of euler angles. This does not make sense.")
                exit(1)
            (euler_angle_0,euler_angle_1,euler_angle_2) = condor.utils.diffraction.random_euler_angles()
        elif self._alignment == "euler_angles":
            # Many orientations (lists of euler angles)
            if isinstance(self._euler_angle_0,list):
                euler_angle_0 = self._euler_angle_0[self._i]
                euler_angle_1 = self._euler_angle_1[self._i]
                euler_angle_2 = self._euler_angle_2[self._i]
            # One orientation (euler angles are scalars)
            else:
                euler_angle_0 = self._euler_angle_0
                euler_angle_1 = self._euler_angle_1
                euler_angle_2 = self._euler_angle_2
        elif self._alignment == "random_euler_angle_0":
            if self._euler_angle_0 is not None:
                log(condor.CONDOR_logger.error,"Conflict of arguments: Specified random_euler_angle_0 alignment and also specified a specific euler_angle_0 = %f. This does not make sense." % self._euler_angle_0)
                exit(1)
            euler_angle_0 = numpy.random.uniform(0,2*numpy.pi)
            euler_angle_1 = self._euler_angle_1 if self._euler_angle_1 is not None else 0.
            euler_angle_2 = self._euler_angle_2 if self._euler_angle_2 is not None else 0.
        return euler_angle_0,euler_angle_1,euler_angle_2

    def set_position_variation(self, position_variation, position_spread, position_variation_n):
        self._position_variation = Variation(position_variation,position_spread,position_variation_n,number_of_dimensions=3,name="particle position")

    def get_position_variation(self):
        A = {
            "position_variation":        self._position_variation.get_mode(),
            "position_variation_spread": self._position_variation.get_spread(),
            "position_variation_n":      self._position_variation.n
        }
        return A
        
    def _get_next_position(self):
        return self._position_variation.get(self.position_mean)

class AbstractContinuousParticle(AbstractParticle):
    def __init__(self,
                 diameter, diameter_variation = None, diameter_spread = None, diameter_variation_n = None,
                 alignment = None, euler_angle_0 = None, euler_angle_1 = None, euler_angle_2 = None,
                 concentration = 1.,
                 position = None,  position_variation = None, position_spread = None, position_variation_n = None,
                 material_type = None, massdensity = None, **atomic_composition):
        
        # Initialise base class
        AbstractParticle.__init__(self,
                                  alignment=alignment, euler_angle_0=euler_angle_0, euler_angle_1=euler_angle_1, euler_angle_2=euler_angle_2,
                                  concentration=concentration,
                                  position=position, position_variation=position_variation, position_spread=position_spread, position_variation_n=position_variation_n)
        # Diameter
        self.set_diameter_variation(diameter_variation=diameter_variation, diameter_spread=diameter_spread, diameter_variation_n=diameter_variation_n)
        self.diameter_mean = diameter
        # Material
        self.set_material(material_type=material_type, massdensity=massdensity, **atomic_composition)

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

    def set_diameter_variation(self,diameter_variation=None,diameter_spread=None,diameter_variation_n=None,**kwargs):
        self._diameter_variation = Variation(diameter_variation,diameter_spread,diameter_variation_n,name="sample diameter")       

    def _get_next_diameter(self):
        d = self._diameter_variation.get(self.diameter_mean)
        # Non-random diameter
        if self._diameter_variation._mode in [None,"range"]:
            if d <= 0:
                log(condor.CONDOR_logger.error,"Sample diameter smaller-equals zero. Change your configuration.")
            else:
                return d
        # Random diameter
        else:
            if d <= 0.:
                log(condor.CONDOR_logger.warning,"Sample diameter smaller-equals zero. Try again.")
                return self._get_next_diameter()
            else:
                return d

    def set_material(self, material_type = None, massdensity = None, **atomic_composition):
        self.material = Material(material_type=material_type, massdensity=massdensity, **atomic_composition)


