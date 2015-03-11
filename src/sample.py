# ----------------------------------------------------------------------------------------------------- 
# CONDOR 
# Simulator for diffractive single-particle imaging experiments with X-ray lasers
# http://xfel.icm.uu.se/condor/
# ----------------------------------------------------------------------------------------------------- 
# Copyright 2014 Max Hantke, Filipe R.N.C. Maia, Tomas Ekeberg
# Condor is distributed under the terms of the GNU General Public License
# ----------------------------------------------------------------------------------------------------- 
# General note:
#  All variables are in SI units by default. Exceptions explicit by variable name.
# ----------------------------------------------------------------------------------------------------- 

import sys,numpy
import scipy.stats
import logging
logger = logging.getLogger("Condor")
if "utils" not in sys.path: sys.path.append("utils")
from variation import Variation
import sample_species

class Sample:
    def __init__(self,**kwargs):
        self.number_of_images = kwargs["number_of_images"]
        self.number_of_particles_mean = kwargs["number_of_particles"]
        self.set_number_of_particles_variation(kwargs["number_of_particles_variation"],kwargs.get("number_of_particles_spread",None),kwargs.get("number_of_particles_variation_n",None))
        self.species = []

    def set_number_of_particles_variation(self,mode,spread,variation_n):
        self._number_of_particles_variation = Variation(mode,spread,variation_n)

    def _get_next_number_of_particles(self):
        N = self._number_of_particles_variation.get(self.number_of_particles_mean)
        # Non-random
        if self._number_of_particles_variation._mode in [None,"range"]:
            if N <= 0:
                logger.error("Sample number of particles smaller-equals zero. Change your configuration.")
                exit(0)
            else:
                return N
        # Random
        else:
            if N <= 0.:
                logger.warning("Sample number of particles smaller-equals zero. Try again.")
                return self.get_next_number_of_particles()
            else:
                return N       

    def add_species(self,**kwargs):
        t = kwargs["sample_type"]
        if t == "uniform_sphere":
            S = sample_species.SampleSpeciesSphere(**kwargs)
        elif t == "uniform_spheroid":
            S = sample_species.SampleSpeciesSpheroid(**kwargs)
        elif t == "map3d":
            S = sample_species.SampleSpeciesMap(**kwargs)
        else:
            logger.error("Sample species for sample_type=%s is not implemented.",t)
            exit(0)
        self.species.append(S)
        
    def _next_particles(self):
        N = self._get_next_number_of_particles()
        if N > 1:
            if len(self.species) == 1:
                self.particles = self.species * N
            else:
                i_s = range(len(self.species))
                c_s = [s.concentration for s in self.species]
                dist = scipy.stats.rv_discrete(name='species distribution', values=(i_s, c_s))
                self.particles = self.species[dist.rvs(size=self.number_of_particles)]
        else:
            self.particles = self.species

    def get_next(self):
        self._next_particles()
        O = {}
        O["sample_types"] = []
        O["particles"] = []
        for p in self.particles:
            O["particles"].append(p.get_next())
        O["number_of_particles"] = len(self.particles)
        return O

