#!/usr/bin/env python

# Installation of propagator
import sys, os, fileinput
import constants_data.fetchsf as sf
from distutils.core import setup

# Check if python_tools are installed
import python_tools.imgtools, python_tools.gentools, python_tools.cxitools
try:
    print "Checking if python_tools are installed..."
    import python_tools.imgtools, python_tools.gentools, python_tools.cxitools
    print "Necessary python tools are installed."
except:
    print "ERROR: Cannot import python_tools. Please install Max' python_tools and add them to your PYTHONPATH before you proceed."
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

print 'Wrapping NFFT...'
pdir = os.path.dirname(os.path.realpath(__file__))
os.chdir("%s/propagator/utils/nfft" % pdir)
os.system("python setup.py build")
os.chdir(pdir)

setup(name='propagator',
      description='Python tools for image analysis',
      version='0.0',
      author='Max Felix Hantke',
      author_email='maxhantke@gmail.com',
      url='github.com/mhantke/propagator',
      packages=['propagator','propagator.utils','propagator.utils.nfft'],
      package_data={'propagator':['data/*'],'propagator.utils.nfft':['nfft_c.so']},
     )

import propagator

print 'Done.'
