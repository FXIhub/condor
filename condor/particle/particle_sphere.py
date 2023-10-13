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

from .particle_abstract import AbstractContinuousParticle

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
                 number = 1., arrival = "synchronised",
                 position = None, position_variation = None, position_spread = None, position_variation_n = None,
                 material_type = None, massdensity = None, atomic_composition = None, electron_density = None): 

        # Initialise base class
        AbstractContinuousParticle.__init__(self,
                                            diameter=diameter, diameter_variation=diameter_variation, diameter_spread=diameter_spread, diameter_variation_n=diameter_variation_n,
                                            rotation_values=None, rotation_formalism=None, rotation_mode="extrinsic",
                                            number=number, arrival=arrival,
                                            position=position, position_variation=position_variation, position_spread=position_spread, position_variation_n=position_variation_n,
                                            material_type=material_type, massdensity=massdensity, atomic_composition=atomic_composition, electron_density=electron_density)
        
    def get_next(self):
        """
        Iterate the parameters and return them as a dictionary
        """
        O = AbstractContinuousParticle.get_next(self)
        O["particle_model"] = "sphere"
        return O

    def get_conf(self):
        """
        Get configuration in form of a dictionary. Another identically configured ParticleMap instance can be initialised by:

        .. code-block:: python

          conf = P0.get_conf()                 # P0: already existing ParticleSphere instance
          P1 = condor.ParticleSpheroid(**conf) # P1: new ParticleSphere instance with the same configuration as P0  
        """
        return AbstractContinuousParticle.get_conf(self)

    def get_dn(self, photon_wavelength):
        if self.materials is None:
            dn = 0.
        else:
            dn = numpy.array([m.get_dn(photon_wavelength) for m in self.materials]).sum()
        return dn
