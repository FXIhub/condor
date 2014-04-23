#!/usr/bin/env python

# Installation of Penguin
import sys, os, fileinput
import constants_data.fetchsf as sf
from distutils.core import setup

# Check if necessary packages are installed
try:
    print "Checking if numpy is installed..."
    import numpy
    print "Necessary package numpy is installed."
except:
    print "ERROR: Cannot import numpy. Please install it and try to install Penguin again."
    print "Installation of Penguin failed."
    quit(0)

try:
    print "Checking if scipy is installed..."
    import scipy
    print "Necessary package scipy is installed."
except:
    print "ERROR: Cannot import scipy. Please install it and try to install Penguin again."
    print "Installation of Penguin failed."
    quit(0)

try:
    print "Checking if h5py is installed..."
    import h5py
    print "Necessary package h5py is installed."
except:
    print "ERROR: Cannot import h5py. Please install it and try to install Penguin again."
    print "Installation of Penguin failed."
    quit(0)

try:
    print "Checking if python_tools are installed..."
    import python_tools.imgtools, python_tools.gentools, python_tools.cxitools
    print "Necessary python tools are installed."
except:
    print "ERROR: Cannot import python_tools. Please install Max' python_tools and add them to your PYTHONPATH before you proceed. You can clone python_tools from git@bitbucket.org:maxhantke/python_tools.git. Execute 'cd ~/target/directory; git clone git@bitbucket.org:maxhantke/python_tools.git; cd python_tools; python setup.py install'."
    print "Installation of Penguin failed."
    quit(0)


# Create dir for data
print 'Clean up data directory'
os.system("mkdir -p ./penguin/data/")

# Copy default configuration file there
print 'Link default configuration file'
os.system("ln -s ./conf/default.conf ./penguin/data/default.conf")

# Scattering factors from the Henke tables and atomic masses 
print 'Loading scattering constants...'
sf.generate_datafile("constants_data/sf","./penguin/data")
print 'Done.'

print 'Wrapping NFFT...'
pdir = os.path.dirname(os.path.realpath(__file__))
os.chdir("%s/penguin/utils/nfft" % pdir)
os.system("python setup.py build")
os.chdir(pdir)

print 'Wrapping ICOSAHEDRON'
os.chdir("%s/penguin/utils/icosahedron" % pdir)
os.system("python setup.py build")
os.chdir(pdir)

setup(name='penguin',
      description='Simulator for diffractive single-particle imaging experiments with X-ray lasers',
      version='1.0',
      author='Max Felix Hantke, Filipe R.N.C. Maia, Tomas Ekeberg',
      author_email='hantke@xray.bmc.uu.se',
      url='http://xfel.icm.uu.se/penguin/',
      packages=['penguin','penguin.utils','penguin.utils.nfft','penguin.utils.icosahedron'],
      package_data={'penguin':['data/*'],'penguin.utils.nfft':['nfft.so'],'penguin.utils.icosahedron':['icosahedron.so']},
     )

# test import
import penguin

print 'Penguin installation was successful.'
