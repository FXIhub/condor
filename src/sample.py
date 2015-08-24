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

import sys,os
import numpy
import scipy.stats

import logging
logger = logging.getLogger(__name__)

import condor.utils.log
from condor.utils.log import log_and_raise_error,log_warning,log_info,log_debug
import condor.utils.tools
from condor.utils.config import load_config
from condor.utils.variation import Variation

def load_sample(conf):
    """
    Create new Sample instance and load parameters from a Condor configuration file or dictionary.
    
    Kwargs:
       :conf(str): Condor configuration file or dictionary (default = None)
    """
    C = condor.utils.config.load_config({"sample": load_config(conf)["sample"]}, {"sample": load_config(condor.CONDOR_default_conf)["sample"]})
    sample = Sample(**C["sample"])
    return sample
        
class Sample:
    """
    Class for sample
    """
    def __init__(self, number_of_particles=1, number_of_particles_variation=None, number_of_particles_spread=None, number_of_particles_variation_n=None, particle_pick="random"):
        """
        Initialisation of a Sample instance. Call load_sample(conf) for initialisation from file 

        Kwargs:
           :number_of_particles(int): Number of particles present at the same time in the interaction volume (default = 1)
           :number_of_particles_variation(str): Variation of the number of particles, either \"poisson\", \"uniform\", \"range\" or None (default = None)
           :number_of_particles_spread(int): Statistical spread of the number of particles (default = None)
           :numbrt_of_particles_variation_n(int): This parameter is only effective if number_of_particles_variation_n=\"range\". In that case this parameter determines the number of samples within the interval number_of_particles +/- number_of_particles_spread/2 (default = None)
           :particle_pick(str): This parameter only takes effect if more than one particle is defined. If the number of particles in the interaction volume exceeds one particles are drawn according to the rule defined by this argument, either \"sequential\" or \"random\" (default = \"random\")
        """
        self.number_of_particles_mean = number_of_particles
        self.particle_pick = particle_pick
        self.set_number_of_particles_variation(number_of_particles_variation, number_of_particles_spread, number_of_particles_variation_n)
        self._particle_models = []
        self._particle_models_names = []

    def get_conf(self):
        """
        Rerieve configuration of this sample instance in form of a dictionary
        """
        conf = {}
        conf["number_of_particles"]             = self.number_of_particles_mean
        conf["number_of_particles_variation"]   = self._number_of_particles_variation.get_mode()
        conf["number_of_particles_spread"]      = self._number_of_particles_variation.get_spread()
        conf["number_of_particles_variation_n"] = self._number_of_particles_variation.n
        conf["particle_pick"]                   = self.particle_pick
        return conf
        
    def get_next(self):
        """
        Draw next
        """
        self._next_particles()
        O = {}
        #O["particle_types"] = []
        O["particles"] = []
        for p in self._particles:
            O["particles"].append(p.get_next())
        O["number_of_particles"] = len(self._particles)
        return O

    def append_particle(self, particle_model, name):
        """
        Add a particle model to the sample
        
        Args:
           :particle_model: Particle model instance
           :name(str): Name of the new particle
        """
        if not isinstance(particle_model, condor.particle.ParticleSphere) and \
           not isinstance(particle_model, condor.particle.ParticleSpheroid) and \
           not isinstance(particle_model, condor.particle.ParticleMap) and \
           not isinstance(particle_model, condor.particle.ParticleMolecule):
            log_and_raise_error(logger, "The argument is not a valid particle model instance.")
        elif name in self._particle_models_names:
            log_and_raise_error(logger, "The given name for the particle model already exists. Please choose andother name and try again.")
        else:
            self._particle_models.append(particle_model)
            self._particle_models_names.append(name)

    def remove_particle(self, name):
        """
        Remove particle model from sample

        Args:
           :name(str): Name of the particle model that shall be removed
        """
        if name not in self._particle_models_names:
            log_and_raise_error(logger, "The given name for the particle model does not exist. Please give an existing particle model name and try again.")
        else:
            i = self._particle_models_names.index(name)
            self._particle_models_names.pop(i)
            self._particle_models.pop(i)

    def remove_all_particles(self):
        """
        Remove all particle models from sample
        """
        self._particle_models = []
        self._particle_models_names = []

    def get_particles(self):
        """
        Get all particle models and names in form of a dictionary
        """
        pm = {}
        for n,p in zip(self._particle_models_names, self._particle_models):
            pm[n] = p
        return pm
            
    def set_number_of_particles_variation(self, mode, spread, variation_n):
        """
        Set statistical variation model for the number of particles

        Args:
           :mode(str): Either None, \"poisson\", \"normal\", \"uniform\" or \"range\"
           :spread(float): Width of the distribution of number of particles. Only matters if mode is \"normal\", \"uniform\" or \"range\". For \'mode\'=\"normal\" \'spread\' means standard deviation, for \'mode\'=\"uniform\" and '\mode'\=\"range\" \'spread\' determines the full width of the distribution
           :variation_n(int): If mode is \"range\" this argument specifies the number of samples within the specified range
        """
        self._number_of_particles_variation = Variation(mode, spread, variation_n)

    def _get_next_number_of_particles(self):
        N = self._number_of_particles_variation.get(self.number_of_particles_mean)
        # Non-random
        if self._number_of_particles_variation._mode in [None,"range"]:
            if N <= 0:
                log_and_raise_error(logger, "Sample number of particles smaller-equals zero. Change your configuration.")
                sys.exit(0)
            else:
                return N
        # Random
        else:
            if N <= 0.:
                log_warning(logger, "Sample number of particles smaller-equals zero. Trying again.")
                return self._get_next_number_of_particles()
            else:
                return N       
        
    def _next_particles(self):
        N = self._get_next_number_of_particles()
        if self.particle_pick == "sequential":
            self._particles = [self._particle_models[i%len(self._particle_models)] for i in range(N)]
        elif self.particle_pick == "random":
            i_s = range(len(self._particle_models))
            c_s = numpy.array([p.concentration for p in self._particle_models])
            c_s = c_s / c_s.sum()
            dist = scipy.stats.rv_discrete(name='model distribution', values=(i_s, c_s))
            self._particles = [self._particle_models[i] for i in dist.rvs(size=N)]
        else:
            log_and_raise_error(logger, "particle_pick=%s is not a valid configuration." % self.particle_variation)
            sys.exit(0)

        
