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

import condor as p
import configobj
import matplotlib,os,numpy,sys
from python_tools import gentools

#pdir = os.path.abspath(os.path.dirname(__file__))+"/../"
#odir = os.path.abspath(os.path.dirname(__file__))+"/ensemble_out/"
#os.system("mkdir -p %s" % odir)

numpy.random.seed(0)

if len(sys.argv) < 2:
    print "ERROR: Please specify a configuration file."
    print "Usage: ./ensemble.py ensemble.conf"
    exit(1)

I = p.Input(sys.argv[1])
O = p.Output(I)
O.write("ensemble.cxi",output=["intensity_pattern","fourier_space_image","real_space_image","sample_diameter","bitmask_image","mask_image","intensity_pattern_center","intensity"])
