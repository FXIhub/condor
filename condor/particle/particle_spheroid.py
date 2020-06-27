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
import numpy

import logging
logger = logging.getLogger(__name__)

import condor
import condor.utils.log
from condor.utils.log import log_and_raise_error,log_warning,log_info,log_debug

from .particle_abstract import AbstractContinuousParticle

from condor.utils.variation import Variation


class ParticleSpheroid(AbstractContinuousParticle):
    """
    Class for a particle model

    *Model:* Uniformly filled spheroid particle (continuum approximation)

    :math:`a`: radius (*semi-diameter*) perpendicular to the rotation axis of the ellipsoid
    :math:`c`: radius (*semi-diameter*) along the rotation axis of the ellipsoid
    
    Before applying rotations the rotation axis is parallel to the the *y*-axis

    Args:
      :diameter (float): Sphere diameter

    Kwargs:
      :diameter_variation (str): See :meth:`condor.particle.particle_abstract.AbstractContinuousParticle.set_diameter_variation` (default ``None``)

      :diameter_spread (float): See :meth:`condor.particle.particle_abstract.AbstractContinuousParticle.set_diameter_variation` (default ``None``)

      :diameter_variation_n (int): See :meth:`condor.particle.particle_abstract.AbstractContinuousParticle.set_diameter_variation` (default ``None``)

      :flattening (float): (Mean) value of :math:`a/c` (default ``0.75``)

      :flattening_variation (str): See :meth:`condor.particle.particle_spheroid.set_flattening_variation` (default ``None``)

      :flattening_spread (float): See :meth:`condor.particle.particle_spheroid.set_flattening_variation` (default ``None``)
    
      :flattening_variation_n (int): See :meth:`condor.particle.particle_spheroid.set_flattening_variation` (default ``None``)

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
                 diameter,
                 diameter_variation = None, diameter_spread = None, diameter_variation_n = None,
                 flattening = 0.75, flattening_variation = None, flattening_spread = None, flattening_variation_n = None,
                 rotation_values = None, rotation_formalism = None, rotation_mode = "extrinsic",
                 number = 1., arrival = "synchronised",
                 position = None, position_variation = None, position_spread = None, position_variation_n = None,
                 material_type = 'water', massdensity = None, atomic_composition = None, electron_density = None):

        # Initialise base class
        AbstractContinuousParticle.__init__(self,
                                            diameter=diameter, diameter_variation=diameter_variation, diameter_spread=diameter_spread, diameter_variation_n=diameter_variation_n,
                                            rotation_values=rotation_values, rotation_formalism=rotation_formalism, rotation_mode=rotation_mode,
                                            number=number, arrival=arrival,
                                            position=position, position_variation=position_variation, position_spread=position_spread, position_variation_n=position_variation_n,
                                            material_type=material_type, massdensity=massdensity, atomic_composition=atomic_composition, electron_density=electron_density)
        self.flattening_mean = flattening
        self.set_flattening_variation(flattening_variation=flattening_variation, flattening_spread=flattening_spread, flattening_variation_n=flattening_variation_n)

    def get_conf(self):
        """
        Get configuration in form of a dictionary. Another identically configured ParticleMap instance can be initialised by:

        .. code-block:: python

          conf = P0.get_conf()                 # P0: already existing ParticleSpheroid instance
          P1 = condor.ParticleSpheroid(**conf) # P1: new ParticleSpheroid instance with the same configuration as P0  
        """
        conf = {}
        conf.update(AbstractContinuousParticle.get_conf(self))
        conf["flattening"] = self.flattening_mean
        fvar = self._flattening_variation.get_conf()
        conf["flattening_variation"] = fvar["mode"]
        conf["flattening_spread"] = fvar["spread"]
        conf["flattening_variation_n"] = fvar["n"]
        return conf
        
    def get_next(self):
        """
        Iterate the parameters and return them as a dictionary
        """
        O = AbstractContinuousParticle.get_next(self)
        O["particle_model"] = "spheroid"
        O["flattening"] = self._get_next_flattening()
        return O
        
    def set_flattening_variation(self, flattening_variation, flattening_spread, flattening_variation_n):
        """
        Set the variation scheme of the flattening parameter
        
        Args:
          :flattening_variation (str): Variation of the particle flattening

            *Choose one of the following options:*
              
              - ``None`` - No variation

              - ``\'normal\'`` - Normal (*Gaussian*) variation

              - ``\'uniform\'`` - Uniformly distributed flattenings

              - ``\'range\'`` - Equidistant sequence of particle-flattening samples within the spread limits. ``flattening_variation_n`` defines the number of samples within the range

          :flattening_spread (float): Statistical spread of the parameter

          :flattening_variation_n (int): Number of particle-flattening samples within the specified range

            .. note:: The argument ``flattening_variation_n`` takes effect only if ``flattening_variation=\'range\'``
        """
        self._flattening_variation = Variation(flattening_variation, flattening_spread, flattening_variation_n)       

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

    def get_dn(self, photon_wavelength):
        if self.materials is None:
            dn = 0.
        else:
            dn = numpy.array([m.get_dn(photon_wavelength) for m in self.materials]).sum()
        return dn
