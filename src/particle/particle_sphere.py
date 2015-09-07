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

from particle_abstract import AbstractContinuousParticle

class ParticleSphere(AbstractContinuousParticle):
    """
    Class for a particle model

    *Model:* Uniformly filled spherical particle (continuum approximation)

    Args:
      :diameter (float): Sphere diameter

    Kwargs:
      :diameter_variation (str): See :meth:`condor.particle.particle_abstract.AbstractContinuousParticle.set_diameter_variation` (default ``None``)

      :diameter_spread (float): See :meth:`condor.particle.particle_abstract.AbstractContinuousParticle.set_diameter_variation` (default ``None``)

      :diameter_variation_n (int): See :meth:`condor.particle.particle_abstract.AbstractContinuousParticle.set_diameter_variation` (default ``None``)

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
                 diameter, diameter_variation = None, diameter_spread = None, diameter_variation_n = None,
                 concentration = 1.,
                 position = None, position_variation = None, position_spread = None, position_variation_n = None,
                 material_type = None, massdensity = None, atomic_composition = None): 

        # Initialise base class
        AbstractContinuousParticle.__init__(self,
                                            diameter=diameter, diameter_variation=diameter_variation, diameter_spread=diameter_spread, diameter_variation_n=diameter_variation_n,
                                            rotation_values=None, rotation_formalism=None, rotation_mode="extrinsic",
                                            concentration=concentration,
                                            position=position, position_variation=position_variation, position_spread=position_spread, position_variation_n=position_variation_n,
                                            material_type=material_type, massdensity=massdensity, atomic_composition=atomic_composition)
        
    def get_next(self):
        """
        Iterate the parameters and return them as a dictionary
        """
        O = AbstractContinuousParticle.get_next(self)
        O["particle_model"] = "sphere"
        return O

    def get_conf(sef):
        """
        Get configuration in form of a dictionary. Another identically configured ParticleMap instance can be initialised by:

        .. code-block:: python

          conf = P0.get_conf()                 # P0: already existing ParticleSphere instance
          P1 = condor.ParticleSpheroid(**conf) # P1: new ParticleSphere instance with the same configuration as P0  
        """
        return AbstractContinuousParticle.get_conf(self)
