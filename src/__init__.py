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

from .experiment import Experiment
from .source import Source, load_source
from .sample import Sample, load_sample
from .particle import ParticleSphere, ParticleSpheroid, ParticleMap, ParticleMolecule, load_particle, load_particles
from .detector import Detector, load_detector

# Set global variables
def _init():
    # Installation directory of Condor
    import os
    global CONDOR_directory
    CONDOR_directory = os.path.dirname(os.path.realpath(__file__))

    # Log to stdout
    import logging, sys
    h = logging.StreamHandler(sys.stdout)
    logging.getLogger('condor').addHandler(h)
    
    # Load data into global variables
    import _data
    data_dir = CONDOR_directory + "/data"
    global CONDOR_atomic_scattering_factors
    CONDOR_atomic_scattering_factors = _data.load_atomic_scattering_factors(data_dir)
    global CONDOR_atomic_masses
    CONDOR_atomic_masses             = _data.load_atomic_masses(data_dir)
    global CONDOR_atomic_numbers
    CONDOR_atomic_numbers            = _data.load_atomic_numbers(data_dir)
    global CONDOR_atomic_compositions
    CONDOR_atomic_compositions       = _data.load_atomic_compositions()
    global CONDOR_mass_densities
    CONDOR_mass_densities            = _data.load_mass_densities()

    # Default configuration files
    global CONDOR_default_pdb
    CONDOR_default_pdb = data_dir + "/DNA.pdb"
    global CONDOR_default_conf
    CONDOR_default_conf = data_dir + "/condor.conf"


_init()

# Delete references to modules for clarity on global scope
#del experiment
#del source
#del sample
#del particle
#del detector
#del utils
