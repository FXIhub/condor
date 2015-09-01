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
    """
    Class for a particle model

    *Model:* Uniformly filled spheroid particle (continuum approximation)

    :math:`a`: radius (*semi-diameter*) perpendicular to the rotation axis of the ellipsoid
    :math:`c`: radius (*semi-diameter*) along the rotation axis of the ellipsoid
    
    Before applying rotations the rotation axis is parallel to the the *y*-axis

    **Arguments:**

      :diameter (float): Sphere diameter

    **Keyword arguments:**
    
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

      :concentration (array): See :class:`condor.particle.particle_abstract.AbstractParticle` (default ``None``)

      :position (array): See :class:`condor.particle.particle_abstract.AbstractParticle` (default ``None``)

      :position_variation (str): See :meth:`condor.particle.particle_abstract.AbstractParticle.set_position_variation` (default ``None``)

      :position_spread (float): See :meth:`condor.particle.particle_abstract.AbstractParticle.set_position_variation` (default ``None``)

      :position_variation_n (int): See :meth:`condor.particle.particle_abstract.AbstractParticle.set_position_variation` (default ``None``)

      :material_type (str): See :meth:`condor.particle.particle_abstract.AbstractContinuousParticle.set_material` (default ``\'water\'``)

      :massdensity (float): See :meth:`condor.particle.particle_abstract.AbstractContinuousParticle.set_material` (default ``None``)

      :atomic_composition (dict): See :meth:`condor.particle.particle_abstract.AbstractContinuousParticle.set_material` (default ``None``)
    """
    def __init__(self,
                 diameter,
                 diameter_variation = None, diameter_spread = None, diameter_variation_n = None,
                 flattening = 0.75, flattening_variation = None, flattening_spread = None, flattening_variation_n = None,
                 rotation_values = None, rotation_formalism = None, rotation_mode = "extrinsic",
                 concentration = 1.,
                 position = None, position_variation = None, position_spread = None, position_variation_n = None,
                 material_type = 'water', massdensity = None, atomic_composition = None):

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
