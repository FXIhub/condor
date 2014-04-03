#!/usr/bin/env python

# Installation of propagator
import sys, os, fileinput
import constants_data.fetchsf as sf
from distutils.core import setup

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


# Check if python_tools are installed
# try:
#     print "Checking if python_tools are installed..."
#     import python_tools.imgtools, python_tools.gentools, python_tools.cxitools
#     print "Necessary python tools are installed."
# except:
#     print "ERROR: Cannot import python_tools. Please install Max' python_tools and add them to your PYTHONPATH before you proceed. You can clone python_tools from git@bitbucket.org:maxhantke/python_tools.git. Execute 'cd ~/target/directory; git clone git@bitbucket.org:maxhantke/python_tools.git; cd python_tools; python setup.py install'."
#     print "Installation of propagator failed."
#     quit(0)

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

print 'Wrapping ICOSAHEDRON'
os.chdir("%s/propagator/utils/icosahedron" % pdir)
os.system("python setup.py build")
os.chdir(pdir)

setup(name='propagator',
      description='Python tools for image analysis',
      version='0.0',
      author='Max Felix Hantke',
      author_email='maxhantke@gmail.com',
      url='github.com/mhantke/propagator',
      packages=['propagator','propagator.utils','propagator.utils.nfft','propagator.utils.icosahedron'],
      package_data={'propagator':['data/*'],'propagator.utils.nfft':['nfft.so'],'propagator.utils.icosahedron':['icosahedron.so']},
     )

# test import
import propagator

print 'Propagator installation successful.'
