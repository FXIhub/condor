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

import logging
logger = logging.getLogger("Condor")
import utils.log
from utils.log import log 

from particle_abstract import AbstractParticleModel

class ParticleModelMolecule(AbstractParticleModel):
    def __init__(self,**kwargs):
        try:
            import spsim
        except:
            log(logger.error,"Cannot import spsim module. This module is necessary to simulate diffraction for particle model \"molecule\". Please install spsim from https://github.com/FilipeMaia/spsim abnd try again.")
            return
        # Check for valid set of keyword arguments
        self.req_keys = [["pdb_filename",["atomic_position","atomic_number"]]]
        self.opt_keys = []
        # Start initialisation
        AbstractParticleModel.__init__(self,**kwargs)
        self.atomic_position   = None
        self.atomic_number     = None
        self.pdb_filename      = None
        if "pdb_filename" in kwargs:
            if kwargs["pdb_filename"] is not None:
                if os.path.isfile(kwargs["pdb_filename"]):
                    self.pdb_filename = kwargs["pdb_filename"]
                else:
                    log(logger.error,"Cannot initialize particle model molecule. PDB file %s does not exist." % kwargs["pdb_filename"])
                    sys.exit(0)
        if self.pdb_filename is None:
            self.atomic_position = numpy.array(kwargs["atomic_position"])
            self.atomic_number   = numpy.array(kwargs["atomic_number"])

    def get_next(self):
        O = AbstractParticleModel.get_next(self)
        O["pdb_filename"]     = self.pdb_filename
        O["atomic_position"] = self.atomic_position
        O["atomic_number"]   = self.atomic_number
        return O


def get_spsim_conf(D_source,D_particle,D_detector):
    s = []
    s += "# THIS FILE WAS CREATED AUTOMATICALLY BY CONDOR\n"
    s += "# Temporary configuration file for spsim\n"
    s += "verbosity_level = 0;\n"
    s += "number_of_dimensions = 2;\n"
    s += "number_of_patterns = 1;\n"
    s += "input_type = \"pdb\";\n"
    s += "pdb_filename = \"%s\";\n" % D_particle["pdb_filename"]
    s += "detector_distance = %.6e;\n" % D_detector["distance"]
    s += "detector_width = %.6e;\n" % (D_detector["pixel_size"] * D_detector["nx"]) 
    s += "detector_height = %.6e;\n" % (D_detector["pixel_size"] * D_detector["ny"])
    s += "detector_pixel_width = %.6e;\n" % D_detector["pixel_size"]
    s += "detector_pixel_height = %.6e;\n" % D_detector["pixel_size"]
    s += "detector_center_x = %.6e;\n" % (D_detector["pixel_size"] * (D_detector["cx"] - (D_detector["nx"]-1)/2.))
    s += "detector_center_y = %.6e;\n" % (D_detector["pixel_size"] * (D_detector["cy"] - (D_detector["ny"]-1)/2.))
    s += "detector_binning = 1;\n"
    s += "experiment_wavelength = %.6e;\n" % D_source["wavelength"]
    s += "experiment_beam_intensity = %.6e;\n" % D_particle["intensity"]
    return s

