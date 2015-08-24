#!/usr/bin/env python

#--------------------------------------------------------------------------
# CONDOR
# URL: xfel.icm.uu.se/condor
# Copyright 2014 Max Hantke, Filipe R.N.C. Maia, Tomas Ekeberg
# Condor is distributed under the terms of the GNU General Public License
#--------------------------------------------------------------------------

packages = ["numpy","scipy","h5py","spsim"]
installed = []

# Check whether necessary packages are installed
try:
    print "Checking whether numpy is installed..."
    import numpy
    print "OK"
    installed.append("numpy")
except:
    print "FAILED"

try:
    print "Checking whether scipy is installed..."
    import scipy
    print "OK"
    installed.append("scipy")
except:
    print "FAILED"

try:
    print "Checking whether h5py is installed..."
    import h5py
    print "OK"
    installed.append("h5py")
except:
    print "FAILED"

# Check whether optional packages are installed
try:
    print "Checking whether spsim is installed..."
    import spsim
    print "OK"
    installed.append("spsim")
except:
    print "FAILED"


if len(installed) == 4:
    print "SUCCESS: All dependencies met."
elif len(installed) == 3 and "spsim" not in installed:
    print "SUCCESS: All necessary dependencies met."
    print "WARNING: Cannot import spsim. If you want the full functionality of your Condor installation please install spsim. It can be downloaded from https://github.com/FilipeMaia/spsim"
else:
    print "FAILED: The following dependencies are not met:"
    print [p for p in packages if p not in installed]
