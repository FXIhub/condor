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

import config
from material import Material
from utils.variation import Variation


class AbstractParticleSpecies:
    def __init__(self, **kwargs):       
        # Check for valid set of keyword arguments
        self.req_keys += ["particle_species","alignment","euler_angle_0","euler_angle_1","euler_angle_2","position","concentration"]
        self.opt_keys += ["position_variation","position_spread","position_variation_n",
                          "geometry_euler_angle_0","geometry_euler_angle_1","geometry_euler_angle_2"] 
        # Check input
        miss_keys,ill_keys = config.check_input(kwargs.keys(),self.req_keys,self.opt_keys,verbose=True)
        if len(miss_keys) > 0: 
            log(logger.error,"Cannot initialize %s because of missing keyword arguments." % self.__class__.__name__)
            exit(1)
        if len(ill_keys) > 0:
            log(logger.error,"Cannot initialize %s instance because of illegal keyword arguments." % self.__class__.__name__)
            exit(1)

        # Start initialisation
        self.set_alignment(alignment=kwargs["alignment"],euler_angle_0=kwargs["euler_angle_0"],euler_angle_1=kwargs["euler_angle_1"],euler_angle_2=kwargs["euler_angle_2"])
        self.set_position_variation(position_variation=kwargs["position_variation"],position_spread=kwargs.get("position_spread",None),position_variation_n=kwargs.get("position_variation_n",None))
        self.position_mean = kwargs["position"]
        self.concentration = kwargs["concentration"]
        self.geometry_euler_angle_0 = kwargs.get("geometry_euler_angle_0",0.)
        self.geometry_euler_angle_1 = kwargs.get("geometry_euler_angle_1",0.)
        self.geometry_euler_angle_2 = kwargs.get("geometry_euler_angle_2",0.)

    def get_next(self):
        O = {}
        O["_class_instance"] = self 
        euler_angle_0,euler_angle_1,euler_angle_2 = self._get_next_orientation()
        O["euler_angle_0"] = euler_angle_0
        O["euler_angle_1"] = euler_angle_1
        O["euler_angle_2"] = euler_angle_2
        O["geometry_euler_angle_0"] = self.geometry_euler_angle_0
        O["geometry_euler_angle_1"] = self.geometry_euler_angle_1
        O["geometry_euler_angle_2"] = self.geometry_euler_angle_2       
        O["position"] = self._get_next_position()
        return O
    
    def set_random_orientation(self):
        self.set_alignment("random")

    def set_alignment(self,alignment=None,euler_angle_0=None,euler_angle_1=None,euler_angle_2=None):
        if alignment not in [None,"random","euler_angles","first_axis","random_euler_angle_0"]:
            log(logger.error,"Invalid argument for sample alignment specified.")
            return
        self._euler_angle_0 = euler_angle_0
        self._euler_angle_1 = euler_angle_1
        self._euler_angle_2 = euler_angle_2
        if alignment is None and (self._euler_angle_0 is not None and self._euler_angle_1 is not None and self._euler_angle_2 is not None):
            self._alignment = "euler_angles"
        else:
            self._alignment = alignment

    def _get_next_orientation(self):
        if self._alignment == "first_axis":
            # Sanity check
            if self._euler_angle_0 is not None or self._euler_angle_1 is not None or self._euler_angle_2 is not None:
                log(logger.error,"Conflict of arguments: Specified first_axis alignment and also specified set of euler angles. This does not make sense.")
                exit(1)
            euler_angle_0 = 0.
            euler_angle_1 = 0.
            euler_angle_2 = 0.
        elif self._alignment == "random":
            # Sanity check
            if self._euler_angle_0 is not None or self._euler_angle_1 is not None or self._euler_angle_2 is not None:
                log(logger.error,"Conflict of arguments: Specified random alignment and also specified set of euler angles. This does not make sense.")
                exit(1)
            (euler_angle_0,euler_angle_1,euler_angle_2) = utils.diffraction.random_euler_angles()
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
                log(logger.error,"Conflict of arguments: Specified random_euler_angle_0 alignment and also specified a specific euler_angle_0 = %f. This does not make sense." % self._euler_angle_0)
                exit(1)
            euler_angle_0 = numpy.random.uniform(0,2*numpy.pi)
            euler_angle_1 = self._euler_angle_1 if self._euler_angle_1 is not None else 0.
            euler_angle_2 = self._euler_angle_2 if self._euler_angle_2 is not None else 0.
        return euler_angle_0,euler_angle_1,euler_angle_2

    def set_position_variation(self,position_variation,position_spread,position_variation_n):
        self._position_variation = Variation(position_variation,position_spread,position_variation_n,number_of_dimensions=3,name="particle position")

    def _get_next_position(self):
        return self._position_variation.get(self.position_mean)

class AbstractContinuousParticleSpecies(AbstractParticleSpecies):
    def __init__(self, **kwargs):
        # Argument check and initialization of inherited class
        self.req_keys = ["diameter"]
        cel_keys = ["c"+k for k in config.DICT_atomic_number.keys()]
        self.opt_keys = ["diameter_variation","diameter_spread","diameter_variation_n","massdensity","material_type"] + cel_keys
        AbstractParticleSpecies.__init__(self,**kwargs)
        # Continue initialization
        # Diameter
        self.set_diameter_variation(diameter_variation=kwargs["diameter_variation"],diameter_spread=kwargs.get("diameter_spread",None),diameter_variation_n=kwargs.get("diameter_variation_n",None))
        self.diameter_mean = kwargs["diameter"]
        # Material
        materialargs = {}
        if 'massdensity' in kwargs:
            materialargs['massdensity'] = kwargs['massdensity']
            for key in kwargs.keys():
                if key in cel_keys: materialargs[key] = kwargs[key]
        elif "material_type" in kwargs:
            materialargs['material_type'] = kwargs['material_type']
        else:
            log(logger.error,"Illegal material configuration for sample species.")
            exit(0)
        self.set_material(**materialargs)

    def get_next(self):
        O = AbstractParticleSpecies.get_next(self)
        O["diameter"] = self._get_next_diameter()
        return O
    
    def set_diameter_variation(self,diameter_variation=None,diameter_spread=None,diameter_variation_n=None,**kwargs):
        self._diameter_variation = Variation(diameter_variation,diameter_spread,diameter_variation_n,name="sample diameter")       

    def _get_next_diameter(self):
        d = self._diameter_variation.get(self.diameter_mean)
        # Non-random diameter
        if self._diameter_variation._mode in [None,"range"]:
            if d <= 0:
                log(logger.error,"Sample diameter smaller-equals zero. Change your configuration.")
            else:
                return d
        # Random diameter
        else:
            if d <= 0.:
                log(logger.warning,"Sample diameter smaller-equals zero. Try again.")
                return self._get_next_diameter()
            else:
                return d

    def set_material(self, **kwargs):
        self.material = Material(**kwargs)
