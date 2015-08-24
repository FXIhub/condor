#!/usr/bin/env python

#--------------------------------------------------------------------------
# CONDOR
# URL: xfel.icm.uu.se/condor
# Copyright 2014 Max Hantke, Filipe R.N.C. Maia, Tomas Ekeberg
# Condor is distributed under the terms of the GNU General Public License
#--------------------------------------------------------------------------

import os,sys
import shutil

import argparse


from distutils.core import setup, Extension

from distutils.sysconfig import get_python_lib

parser = argparse.ArgumentParser(description='Deinstallation script for Condor')
parser.add_argument('-p', '--prefix', metavar='prefix', type=str, 
                    help="prefix for condor installation")
args = parser.parse_args()


if __name__ == "__main__":

    # Delete build directory
    d = os.path.dirname(os.path.realpath(__file__)) + "/build"
    if os.path.exists(d):
        print "Removing build directory %s ..." % d
        shutil.rmtree("%s/" % d)
        print "Done"
    else:
        print "WARNING: Cannot find build directory in %s. Skipping removal of build directory." % d

    # Delete data files
    d = os.path.dirname(os.path.realpath(__file__)) + "/src/data/"
    for f in [d+f for f in os.listdir(d) if f.endswith(".dat")]: 
        print "Removing data file %s ..." % f
        os.remove(f)
        print "Done"
    
    # Check existence of installation directory
    if args.prefix:
        d = get_python_lib(prefix=args.prefix) + "/condor"
    else:
        d = get_python_lib() + "/condor"
    if not os.path.exists(d):
        print "ERROR:\tCannot find or access Condor installation in %s." % d
        print "\tUse the --prefix argument if your installation is not in the standard installation directory."
    else:
        # Delete installation directory
        print "Removing directory %s ..." % d 
        shutil.rmtree("%s/" % d)
        print "Done"

