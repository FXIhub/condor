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
    """
    Class for a particle model

    *Model:* Discrete atomic positions

    Kwargs:
      :pdb_filename (str): See :meth:`set_atoms_from_pdb_file` (default ``None``)

      :atomic_numbers (array): See :meth:`set_atoms_from_arrays` (default ``None``)
    
      :atomic_positions (array): See :meth:`set_atoms_from_arrays` (default ``None``)

      .. note:: The atomic positions have to be specified either by a ``pdb_filename`` or by ``atomic_numbers`` and  ``atomic_positions``.

      :rotation_values (array): See :meth:`condor.particle.particle_abstract.AbstractParticle.set_alignment` (default ``None``)

      :rotation_formalism (str): See :meth:`condor.particle.particle_abstract.AbstractParticle.set_alignment` (default ``None``)

      :rotation_mode (str): See :meth:`condor.particle.particle_abstract.AbstractParticle.set_alignment` (default ``None``)
    
      :concentration (float): See :class:`condor.particle.particle_abstract.AbstractParticle` (default ``None``)

      :position (array): See :class:`condor.particle.particle_abstract.AbstractParticle` (default ``None``)

      :position_variation (str): See :meth:`condor.particle.particle_abstract.AbstractParticle.set_position_variation` (default ``None``)

      :position_spread (float): See :meth:`condor.particle.particle_abstract.AbstractParticle.set_position_variation` (default ``None``)

      :position_variation_n (int): See :meth:`condor.particle.particle_abstract.AbstractParticle.set_position_variation` (default ``None``)
    """
    def __init__(self,
                 pdb_filename = None,
                 atomic_numbers = None, atomic_positions = None,
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
        if pdb_filename is not None and (atomic_numbers is None and atomic_positions is None):
            if os.path.isfile(pdb_filename):
                self.set_atoms_from_pdb_file(pdb_filename)
            else:
                log_and_raise_error(logger, "Cannot initialize particle model molecule. PDB file %s does not exist." % pdb_filename)
                sys.exit(0)

        elif pdb_filename is None and (atomic_numbers is not None and atomic_positions is not None):
            self.set_atoms_from_arrays(atomic_numbers, atomic_positions)
        else:
            log_and_raise_error(logger, "Cannot initialise particle model molecule. The atomic positions have to be specified either by a pdb_filename or by atomic_numbers and atomic_positions.")

    def get_conf(self):
        """
        Get configuration in form of a dictionary. Another identically configured ParticleMolecule instance can be initialised by:

        .. code-block:: python

          conf = P0.get_conf()                 # P0: already existing ParticleMolecule instance
          P1 = condor.ParticleMolecule(**conf) # P1: new ParticleMolcule instance with the same configuration as P0  
        """
        conf = {}
        conf.update(AbstractParticle.get_conf())
        conf["atomic_numbers"]   = self.get_atomic_numbers()
        conf["atomic_positions"] = self.get_atomic_positions()
        return conf
        
    def set_atoms_from_pdb_file(self, pdb_filename):
        """
        Specify atomic positions from a PDB file 

        The PDB file format is described here: `http://www.wwpdb.org/documentation/file-format <http://www.wwpdb.org/documentation/file-format>`

        Args:
          :pdb_filename (str): Location of the PDB file
        """
        import spsim
        mol = spsim.get_Molecule_from_pdb(pdb_filename)
        self._atomic_numbers, self._atomic_positions = spsim.get_atoms_from_molecule(mol)
        spsim.free_mol(mol)
        
    def set_atoms_from_arrays(self, atomic_numbers, atomic_positions):
        r"""
        Specify atomic positions from atomic numbers and atomic positions

        Args:
          :atomic_numbers (array): Integer array of atomic numbers specifies the element species of each atom. Array shape: (:math:`N`,) with :math:`N` denoting the number of atoms.

          :atomic_position (array): Float array of atomic positions [:math:`x`, :math:`y`, :math:`z`] in unit meter. Array shape: (:math:`N`, 3,) with :math:`N` denoting the number of atoms
        """
        N1 = len(atomic_numbers)
        N2 = len(atomic_positions)
        if N1 != N2:
            log_and_raise_error(logger, "Cannot set atoms. atomic_numbers and atomic_positions have to have the same length")
        self._atomic_positions = numpy.array(atomic_positions)
        self._atomic_numbers   = numpy.array(atomic_numbers)

    def get_atomic_numbers(self):
        """
        Return the array of atomic numbers
        """
        return self._atomic_numbers.copy()

    def get_atomic_positions(self):
        """
        Return the array of atomic positions
        """
        return self._atomic_positions.copy()

    def get_atomic_standard_weights(self):
        """
        Return the atomic standard weights in unified atomic mass unit (*u*)
        """
        Z = self.get_atomic_numbers()
        names = [condor.utils.material.atomic_names[z-1] for z in Z]
        M = numpy.array([condor.utils.material.atomic_masses[n] for n in names], dtype=numpy.float64)
        return M
    
    def get_radius_of_gyration(self):
        r"""
        Return the radius of gyration :math:`R_g`

        Atomic structure of :math:`N` atoms with masses :math:`m_i` at the positions :math:`\vec{r}_i`

        :math:`R_g = \fract{ \sqrt{ \sum_{i=0}^N{ \vec{r}_i-\vec{r}_{\text{COM}} } } }{ \sum_{i=0}^N{ m_i }}`
        """
        M = self.get_atomic_masses()
        r = self.get_atomic_positions()
        r_com = self.get_center_of_mass()
        r_g = numpy.sqrt( (M*(r-r_com)**2).sum() / M.sum() )
        return r_g

    def get_center_of_mass(self):
        r"""
        Return the position of the center of mass :math:`\vec{r}_{\text{COM}}`

        Atomic structure of :math:`N` atoms with masses :math:`m_i` at the positions :math:`\vec{r}_i`

        :math:`\vec{r}_{\text{COM}} = \frac{\sum_{i=0}^N{m_i \, \vec{r}_i}}{\sum_{i=0}^N{m_i}}`
        """
        M = self.get_atomic_masses()
        r = self.get_atomic_positions()
        r_com = (r*M).sum() / M.sum()
        return r_com
            
    @property
    def diameter_mean(self):
        """
        Return the two times the radius of gyration as an estimate for the extent (diameter) of the atomic structure
        """
        self._diameter_mean = 2*self.get_radius_of_gyration()
            
    def get_next(self):
        """
        Iterate the parameters and return them as a dictionary
        """
        O = AbstractParticle.get_next(self)
        O["particle_model"]   = "molecule"
        O["atomic_numbers"]   = self.get_atomic_numbers()
        O["atomic_positions"] = self.get_atomic_positions()
        return O

