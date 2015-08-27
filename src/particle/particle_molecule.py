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

import os,sys
import numpy
import tempfile

import logging
logger = logging.getLogger(__name__)

import condor
import condor.utils.log
from condor.utils.log import log_and_raise_error,log_warning,log_info,log_debug

from particle_abstract import AbstractParticle

class ParticleMolecule(AbstractParticle):
    def __init__(self,
                 pdb_filename = None, atomic_positions = None, atomic_numbers = None,
                 rotation_values = None, rotation_formalism = None, rotation_mode = "extrinsic",
                 concentration = 1.,
                 position = None,  position_variation = None, position_spread = None, position_variation_n = None):
        try:
            import spsim
        except:
            log_and_raise_error(logger, "Cannot import spsim module. This module is necessary to simulate diffraction for particle model \"molecule\". Please install spsim from https://github.com/FilipeMaia/spsim abnd try again.")
            return
        # Initialise base class
        AbstractParticle.__init__(self,
                                  rotation_values=rotation_values, rotation_formalism=rotation_formalism, rotation_mode=rotation_mode,                                            
                                  concentration=concentration,
                                  position=position, position_variation=position_variation, position_spread=position_spread, position_variation_n=position_variation_n)
        self._atomic_positions  = None
        self._atomic_numbers    = None
        self._pdb_filename      = None
        self._diameter_mean    = None
        if pdb_filename is not None:
            if os.path.isfile(pdb_filename):
                self.set_atoms_from_pdb(pdb_filename)
            else:
                log_and_raise_error(logger, "Cannot initialize particle model molecule. PDB file %s does not exist." % pdb_filename)
                sys.exit(0)
        else:
            self.set_atoms_from_arrays(atomic_numbers, atomic_positions)

    def get_conf(self):
        conf = {}
        conf.update(AbstractParticle.get_conf())
        nr,pos = self.get_atoms()
        conf["atomic_positions"] = pos
        conf["atomic_numbers"]   = nr
        return conf
        
    def set_atoms_from_pdb(self, pdb_filename):
        import spsim
        mol = spsim.get_Molecule_from_pdb(pdb_filename)
        self._atomic_numbers, self._atomic_positions = spsim.get_atoms_from_molecule(mol)
        spsim.free_mol(mol)
        
    def set_atoms_from_arrays(self, atomic_numbers, atomic_positions):
        self._atomic_positions = numpy.array(atomic_positions)
        self._atomic_numbers   = numpy.array(atomic_numbers)

    def get_atoms(self):
        return [self._atomic_numbers, self._atomic_positions]
        
    def get_radius_of_gyration(self):
        Z,pos = self.get_atoms()
        r = numpy.sqrt( (Z * (pos-pos.mean(axis=0))**2).sum() / float(Z.sum()) )
        return r
        
    def get_atoms(self):
        return self._atomic_numbers, self._atomic_positions
            
    @property
    def diameter_mean(self):
        self._diameter_mean = 2*self.get_radius_of_gyration()
            
    def get_next(self):
        O = AbstractParticle.get_next(self)
        O["particle_model"] = "molecule"
        O["atomic_numbers"], O["atomic_positions"] = self.get_atoms()
        return O

