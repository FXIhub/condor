#!/usr/bin/env python

# Installation of propagator
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
    print "ERROR: Cannot import numpy. Please install it and try to install propagator again."
    print "Installation of propagator failed."
    quit(0)

try:
    print "Checking if scipy is installed..."
    import scipy
    print "Necessary package scipy is installed."
except:
    print "ERROR: Cannot import scipy. Please install it and try to install propagator again."
    print "Installation of propagator failed."
    quit(0)

try:
    print "Checking if h5py is installed..."
    import h5py
    print "Necessary package h5py is installed."
except:
    print "ERROR: Cannot import h5py. Please install it and try to install propagator again."
    print "Installation of propagator failed."
    quit(0)


# Create dir for data
print 'Clean up data directory'
os.system("mkdir -p ./propagator/data/")

# Copy default configuration file there
print 'Copy default configuration file'
os.system("cp ./conf/default.conf ./propagator/data/")

# Scattering factors from the Henke tables and atomic masses 
print 'Loading scattering constants...'
sf.generate_datafile("constants_data/sf","./propagator/data")
print 'Done.'

setup(name='propagator',
      description='Python tools for image analysis',
      version='0.0',
      author='Max Felix Hantke',
      author_email='maxhantke@gmail.com',
      url='github.com/mhantke/propagator',
      #packages=['propagator','propagator.utils','propagator.utils.nfft','propagator.utils.icosahedron'],
      packages=['propagator', 'propagator.utils', "propagator.utils.python_tools"],
      #package_data={'propagator':['data/*'],'propagator.utils.nfft':['nfft.so'],'propagator.utils.icosahedron':['icosahedron.so']},
      package_data={'propagator':['data/*']},
      ext_modules=[Extension("propagator.utils.icosahedron", ["propagator/utils/icosahedron/icosahedronmodule.c"]),
                   Extension("propagator.utils.nfft", sources=["propagator/utils/nfft/nfftmodule.c"], libraries=["nfft3"])],
      include_dirs=numpy.distutils.misc_util.get_numpy_include_dirs()
     )



