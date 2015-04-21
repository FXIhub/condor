#!/usr/bin/env python

#--------------------------------------------------------------------------
# CONDOR
# URL: xfel.icm.uu.se/condor
# Copyright 2014 Max Hantke, Filipe R.N.C. Maia, Tomas Ekeberg
# Condor is distributed under the terms of the GNU General Public License
#--------------------------------------------------------------------------

# Installation of Condor
import sys, os, fileinput
import constants_data.fetchsf as sf
from distutils.core import setup, Extension
import numpy.distutils.misc_util

# Check if necessary packages are installed
try:
    print "Checking if numpy is installed..."
    import numpy
    print "Necessary package numpy is installed."
except:
    print "ERROR: Cannot import numpy. Please install it and try to install Condor again."
    print "Installation of Condor failed."
    quit(0)

try:
    print "Checking if scipy is installed..."
    import scipy
    print "Necessary package scipy is installed."
except:
    print "ERROR: Cannot import scipy. Please install it and try to install Condor again."
    print "Installation of Condor failed."
    quit(0)

try:
    print "Checking if h5py is installed..."
    import h5py
    print "Necessary package h5py is installed."
except:
    print "ERROR: Cannot import h5py. Please install it and try to install Condor again."
    print "Installation of Condor failed."
    quit(0)
    
# Create dir for data
print 'Clean up data directory'
os.system("mkdir -p ./src/data/")

# Copy default configuration file there
print 'Copy default configuration file'
os.system("cp ./default/*.conf ./src/data/")

# Scattering factors from the Henke tables and atomic masses 
print 'Loading scattering constants...'
sf.generate_datafile("constants_data/sf","./src/data")
print 'Done.'

setup(name='condor',
      description='Simulator for diffractive single-particle imaging experiments with X-ray lasers',
      version='1.0',
      author='Max Felix Hantke, Filipe R.N.C. Maia, Tomas Ekeberg',
      author_email='hantke@xray.bmc.uu.se',
      url='http://xfel.icm.uu.se/condor/',
      package_dir={"condor":"src"},
      packages=['condor', 'condor.utils'],
      package_data={'condor':['data/*']},
      scripts=['bin/condor','bin/test_condor'],
      ext_modules=[Extension("condor.utils.icosahedron", sources=["src/utils/icosahedron/icosahedronmodule.c"]),
                   Extension("condor.utils.nfft", sources=["src/utils/nfft/nfftmodule.c"], libraries=["nfft3"])],
      include_dirs=numpy.distutils.misc_util.get_numpy_include_dirs()
     )
