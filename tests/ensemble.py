#!/usr/bin/env python

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

import condor
import configobj
import matplotlib,os,numpy,sys
from python_tools import gentools

#condor.logger.setLevel("DEBUG")

numpy.random.seed(0)

out_dir = os.path.abspath(os.path.dirname(__file__))+"/ensemble_out"
os.system("mkdir -p %s" % out_dir)
out_fn = out_dir + "/ensemble.cxi"

I = condor.Input("ensemble.conf")

O = condor.Output(I)

O.write(out_fn,output="all")#["intensity_pattern","fourier_space_image","real_space_image","sample_diameter","bitmask_image","mask_image","intensity_pattern_center","intensity"])
