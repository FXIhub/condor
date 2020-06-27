# -----------------------------------------------------------------------------------------------------
# CONDOR
# Simulator for diffractive single-particle imaging experiments with X-ray lasers
# http://xfel.icm.uu.se/condor/
# -----------------------------------------------------------------------------------------------------
# Copyright 2016 Max Hantke, Filipe R.N.C. Maia, Tomas Ekeberg
# Condor is distributed under the terms of the BSD 2-Clause License
# -----------------------------------------------------------------------------------------------------
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# 
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# -----------------------------------------------------------------------------------------------------
# General note:
# All variables are in SI units by default. Exceptions explicit by variable name.
# -----------------------------------------------------------------------------------------------------

from __future__ import print_function, absolute_import # Compatibility with python 2 and 3

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
from condor.utils.material import AtomDensityMaterial, ElectronDensityMaterial
from condor.utils.variation import Variation

import condor.utils.diffraction

class AbstractParticle:
    r"""
    Base class for every derived particle class

    Kwargs:
      :rotation_values: See :meth:`condor.particle.particle_abstract.AbstractParticle.set_alignment` (default ``None``)

      :rotation_formalism (str): See :meth:`condor.particle.particle_abstract.AbstractParticle.set_alignment` (default ``None``)

      :rotation_mode (str): See :meth:`condor.particle.particle_abstract.AbstractParticle.set_alignment` (default ``None``)

      :number (float): Expectation value for the number of particles in the interaction volume. (defaukt ``1.``)

      :arrival (str): Arrival of particles at the interaction volume can be either ``'random'`` or ``'synchronised'``. If ``sync`` at every event the number of particles in the interaction volume equals the rounded value of ``number``. If ``'random'`` the number of particles is Poissonian and ``number`` is the expectation value. (default ``'synchronised'``)
    
      :position: (Mean) position vector [*x*, *y*, *z*] of the particle. If set to ``None`` the particle is placed at the origin (default ``None``)

      :position_variation (str): See :meth:`condor.particle.particle_abstract.AbstractParticle.set_position_variation` (default ``None``)

      :position_spread (float): See :meth:`condor.particle.particle_abstract.AbstractParticle.set_position_variation` (default ``None``)

      :position_variation_n (int): See :meth:`condor.particle.particle_abstract.AbstractParticle.set_position_variation` (default ``None``)

    """
    def __init__(self,
                 rotation_values = None, rotation_formalism = None, rotation_mode = "extrinsic",
                 number = 1., arrival = "synchronised",
                 position = None,  position_variation = None, position_spread = None, position_variation_n = None):
        self.set_alignment(rotation_values=rotation_values, rotation_formalism=rotation_formalism, rotation_mode=rotation_mode)
        self.set_position_variation(position_variation=position_variation, position_spread=position_spread, position_variation_n=position_variation_n)
        self.position_mean = position if position is not None else [0., 0., 0.]
        self.number = number
        self.arrival = arrival

    def get_next_number_of_particles(self):
        """
        Iterate the number of partices
        """
        if self.arrival == "random":
            return int(numpy.random.poisson(self.number))
        elif self.arrival == "synchronised":
            return int(numpy.round(self.number))
        else:
            log_and_raise_error(logger, "self.arrival=%s is invalid. Has to be either \'synchronised\' or \'random\'." % self.arrival)
        
    def get_next(self):
        """
        Iterate the parameters of the Particle instance and return them as a dictionary
        """
        O = {}
        O["_class_instance"]      = self
        O["extrinsic_quaternion"] = self._get_next_extrinsic_rotation().get_as_quaternion()
        O["position"]             = self._get_next_position()
        return O

    def get_current_rotation(self):
        """
        Return current orientation of the particle in form of an instance of :class:`condor.utils.rotation.Rotation`
        """
        return self._rotations.get_current_rotation()

    def set_alignment(self, rotation_values, rotation_formalism, rotation_mode):
        """
        Set rotation scheme of the partice

        Args:        
          :rotation_values: Array of rotation parameters. For simulating patterns of many shots this can be also a sequence of rotation parameters. Input ``None`` for no rotation and for random rotation formalisms. For more documentation see :class:`condor.utils.rotation.Rotations` (default ``None``)  

          :rotation_mode (str): If the rotation shall be assigned to the particle choose ``\'extrinsic\'``. Choose ``\'intrinsic\'`` if the coordinate system shall be rotated (default ``\'extrinsic\'``)

        """
        # Check input
        if rotation_mode not in ["extrinsic","intrinsic"]:
            log_and_raise_error(logger, "%s is not a valid rotation mode for alignment." % rotation_mode)
            sys.exit(1)
        self._rotation_mode = rotation_mode
        self._rotations = condor.utils.rotation.Rotations(values=rotation_values, formalism=rotation_formalism)

    def set_position_variation(self, position_variation, position_spread, position_variation_n):
        r"""
        Set position variation scheme

        Args:
          :position_variation (str): Statistical variation of the particle position (default ``None``)

            *Choose one of the following options:*
        
            ====================== ============================================================================================
            ``position_variation`` Type of variation
            ====================== ============================================================================================
            ``None``               No positional variation
            ``'normal'``           Normal (*Gaussian*) variation
            ``'uniform'``          Uniformly distributed positions within spread limits
            ``'range'``            Equidistant sequence of ``position_variation_n`` position samples within ``position_spread``
            ====================== ============================================================================================

          :position_spread (float): Statistical spread of the particle position

          :position_variation_n (int): Number of position samples within the specified range in each dimension

            .. note:: The argument ``position_variation_n`` takes effect only in combination with ``position_variation='range'``
        """
        self._position_variation = Variation(position_variation,position_spread,position_variation_n,number_of_dimensions=3)

    
    def _get_next_extrinsic_rotation(self):
        rotation = self._rotations.get_next_rotation()
        if self._rotation_mode == "intrinsic":
            rotation = copy.deepcopy(rotation)
            rotation.invert()
        return rotation

    def _get_next_position(self):
        return self._position_variation.get(self.position_mean)
    
    def get_conf(self):
        """
        Get configuration in form of a dictionary
        """
        conf = {}
        conf.update(self._get_conf_rotation())
        conf.update(self._get_conf_position_variation())
        conf["number"] = self.number
        conf["arrival"]        = self.arrival
        return conf

    def _get_conf_alignment(self):
        R = self.get_current_rotation()
        A = {
            "rotation_values"    : self._rotations.get_all_values(),
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
    """
    Base class for derived particle classes that make use of the continuum approximation (density instead of discrete atoms)

    Args:
      :diameter (float): (Mean) particle diameter in unit meter
    
    Kwargs:
      :diameter_variation (str): See :meth:`condor.particle.particle_abstract.AbstractContinuousParticle.set_diameter_variation` (default ``None``)

      :diameter_spread (float): See :meth:`condor.particle.particle_abstract.AbstractContinuousParticle.set_diameter_variation` (default ``None``)

      :diameter_variation_n (int): See :meth:`condor.particle.particle_abstract.AbstractContinuousParticle.set_diameter_variation` (default ``None``)

      :rotation_values (array): See :meth:`condor.particle.particle_abstract.AbstractParticle.set_alignment` (default ``None``)

      :rotation_formalism (str): See :meth:`condor.particle.particle_abstract.AbstractParticle.set_alignment` (default ``None``)

      :rotation_mode (str): See :meth:`condor.particle.particle_abstract.AbstractParticle.set_alignment` (default ``None``)

      :number (float): Expectation value for the number of particles in the interaction volume. (defaukt ``1.``)

      :arrival (str): Arrival of particles at the interaction volume can be either ``'random'`` or ``'synchronised'``. If ``sync`` at every event the number of particles in the interaction volume equals the rounded value of ``number``. If ``'random'`` the number of particles is Poissonian and ``number`` is the expectation value. (default ``'synchronised'``)

      :position (array): See :class:`condor.particle.particle_abstract.AbstractParticle` (default ``None``)

      :position_variation (str): See :meth:`condor.particle.particle_abstract.AbstractParticle.set_position_variation` (default ``None``)

      :position_spread (float): See :meth:`condor.particle.particle_abstract.AbstractParticle.set_position_variation` (default ``None``)

      :position_variation_n (int): See :meth:`condor.particle.particle_abstract.AbstractParticle.set_position_variation` (default ``None``)

      :material_type (str): See :meth:`condor.particle.particle_abstract.AbstractContinuousParticle.set_material` (default ``\'water\'``)

      :massdensity (float): See :meth:`condor.particle.particle_abstract.AbstractContinuousParticle.set_material` (default ``None``)

      :atomic_composition (dict): See :meth:`condor.particle.particle_abstract.AbstractContinuousParticle.set_material` (default ``None``)

      :electron_density (float): See :meth:`condor.particle.particle_abstract.AbstractContinuousParticle.set_material` (default ``None``)
    """
    def __init__(self,
                 diameter, diameter_variation = None, diameter_spread = None, diameter_variation_n = None,
                 rotation_values = None, rotation_formalism = None, rotation_mode = "extrinsic",
                 number = 1., arrival = "synchronised",
                 position = None,  position_variation = None, position_spread = None, position_variation_n = None,
                 material_type = 'water', massdensity = None, atomic_composition = None, electron_density = None):
        
        # Initialise base class
        AbstractParticle.__init__(self,
                                  rotation_values=rotation_values, rotation_formalism=rotation_formalism, rotation_mode=rotation_mode,
                                  number=number, arrival=arrival,
                                  position=position, position_variation=position_variation, position_spread=position_spread, position_variation_n=position_variation_n)
        # Diameter
        self.set_diameter_variation(diameter_variation=diameter_variation, diameter_spread=diameter_spread, diameter_variation_n=diameter_variation_n)
        self.diameter_mean = diameter
        # Material
        self.set_material(material_type=material_type, massdensity=massdensity, atomic_composition=atomic_composition, electron_density=electron_density)

    def get_conf(self):
        """
        Get configuration in form of a dictionary
        """
        conf = {}
        conf.update(AbstractParticle.get_conf(self))
        conf["diameter"] = self.diameter_mean
        dvar = self._diameter_variation.get_conf()
        conf["diameter_variation"] = dvar["mode"]
        conf["diameter_spread"] = dvar["spread"]
        conf["diameter_variation_n"] = dvar["n"]
        conf.update(self._get_material_conf())
        return conf
        
    def get_next(self):
        """
        Iterate the parameters of the Particle instance and return them as a dictionary
        """
        O = AbstractParticle.get_next(self)
        O["diameter"] = self._get_next_diameter()
        return O

    def set_diameter_variation(self, diameter_variation, diameter_spread, diameter_variation_n):
        r"""
        Set the variation scheme of the particle diameter
        
        Args:
          :diameter_variation (str): Variation of the particle diameter

            *Choose one of the following options:*

            ====================== ============================================================================================
            ``diameter_variation`` Type of variation
            ====================== ============================================================================================
            ``None``               No diameter variation
            ``'normal'``           Normal (*Gaussian*) variation
            ``'uniform'``          Uniformly distributed diameters within spread limits
            ``'range'``            Equidistant sequence of ``diameter_variation_n`` diameter samples within ``diameter_spread``
            ====================== ============================================================================================

          :diameter_spread (float): Statistical spread

          :diameter_variation_n (int): Number of particle-diameter samples within the specified range

            .. note:: The argument ``diameter_variation_n`` takes effect only if ``diameter_variation='range'``
        """
        self._diameter_variation = Variation(diameter_variation, diameter_spread, diameter_variation_n)       

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

    def set_material(self, material_type, massdensity, atomic_composition, electron_density):
        """
        Initialise and set the AtomDensityMaterial / ElectronDensityMaterial class instance of the particle

        Args:
          :material_type (str): See :class:`condor.utils.material.AtomDensityMaterial`

          :massdensity (float): See :class:`condor.utils.material.AtomDensityMaterial`

          :atomic_composition (dict): See :class:`condor.utils.material.AtomDensityMaterial`

          :electron_density (float): See :class:`condor.utils.material.ElectronDensityMaterial`
        """
        if material_type is None:
            self.materials = None
        else:
            self.materials = []
            if isinstance(material_type, list) or isinstance(massdensity, list) or isinstance(atomic_composition, list) or isinstance(electron_density, list):
                L = max([len(v) for v in [material_type, massdensity, atomic_composition, electron_density] if isinstance(v, list)])
                material_types      = material_type if material_type is not None else [None]*L
                massdensities       = massdensity if massdensity is not None else [None]*L
                atomic_compositions = atomic_composition if atomic_composition is not None else [None]*L
                electron_densities  = electron_density if electron_density is not None else [None]*L
                for material_type_i, massdensity_i, atomic_composition_i, electron_density_i in zip(material_types, massdensities, atomic_compositions, electron_densities):
                    self.add_material(material_type=material_type_i, massdensity=massdensity_i, atomic_composition=atomic_composition_i, electron_density=electron_density_i)
            else:
                self.add_material(material_type=material_type, massdensity=massdensity, atomic_composition=atomic_composition, electron_density=electron_density)

    def add_material(self, material_type, massdensity, atomic_composition, electron_density):
        """
        Initialise and add the AtomDensityMaterial / ElectronDensityMaterial class instance to the particle

        Args:
          :material_type (str): See :class:`condor.utils.material.AtomDensityMaterial`

          :massdensity (float): See :class:`condor.utils.material.AtomDensityMaterial`

          :atomic_composition (dict): See :class:`condor.utils.material.AtomDensityMaterial`

          :electron_density (float): See :class:`condor.utils.material.ElectronDensityMaterial`
        """
        if electron_density is None:
            self.materials.append(AtomDensityMaterial(material_type=material_type, massdensity=massdensity, atomic_composition=atomic_composition))
        else:
            if massdensity is not None or atomic_composition is not None:
                log_and_raise_error(logger, r"An electron density is defined so material_type, massdensity and atomic_composition have to be all 'None'.")
                return
            if material_type != "custom":
                log_and_raise_error(logger, r"An electron density is defined, so material_type must be \'custom\' but is %s." % material_type)
                return
            self.materials.append(ElectronDensityMaterial(electron_density=electron_density))

    def _get_material_conf(self):
        conf = {}
        for m_i in self.materials:
            conf_i = m_i.get_conf()
            if isinstance(m_i, AtomDensityMaterial):
                conf_i["electron_density"] = None
            elif isinstance(m_i, ElectronDensityMaterial):
                conf_i["material_type"] = None
                conf_i["massdensity"] = None
                conf_i["atomic_composition"] = None
            else:
                log_and_raise_error(logger, "Material has the wrong class: %s" % str(m_i))
            for k,v in conf_i.items():
                if k not in conf:
                    conf[k] = []
                conf[k].append(conf_i[k])
        return conf
