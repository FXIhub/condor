#!/usr/bin/env python

#--------------------------------------------------------------------------
# CONDOR
# URL: xfel.icm.uu.se/condor
# Copyright 2014 Max Hantke, Filipe R.N.C. Maia, Tomas Ekeberg
# Condor is distributed under the terms of the GNU General Public License
#--------------------------------------------------------------------------

import os
from distutils.core import setup, Extension

from distutils.sysconfig import get_python_lib

# Delete installed version
d = get_python_lib() + "/condor"
print "Removing directory %s ..." % d 
os.system("rm -rf %s" % d)
print "Done"

# Delete build directory
d = os.path.dirname(os.path.realpath(__file__)) + "/build"
print "Removing build directory %s ..." % d
os.system("rm -rf %s" % d)
print "Done"
