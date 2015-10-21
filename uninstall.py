#!/usr/bin/env python

#--------------------------------------------------------------------------
# CONDOR
# URL: xfel.icm.uu.se/condor
# Copyright 2014 Max Hantke, Filipe R.N.C. Maia, Tomas Ekeberg
# Condor is distributed under the terms of the GNU General Public License
#--------------------------------------------------------------------------

import os,sys
import time
import shutil

import argparse


from distutils.core import setup, Extension

from distutils.sysconfig import get_python_lib

parser = argparse.ArgumentParser(description='Deinstallation script for Condor')
parser.add_argument('-p', '--prefix', metavar='prefix', type=str, 
                    help="prefix for condor installation")
args = parser.parse_args()


if __name__ == "__main__":
   
    # Check existence of installation directory
    if args.prefix:
        plib = get_python_lib(prefix=args.prefix)
    else:
        plib = get_python_lib()
    dirs = [(plib + "/" + d) for d in os.listdir(plib) if d.startswith("condor")]
        
    if len(dirs) == 0:
        try:
            import condor
            dirs = [os.path.dirname(os.path.realpath(condor.__file__))]
        except:
            print "WARNING:\tCannot find or access Condor installation"
            print "\t\tUse the --prefix argument if your installation is not in the standard installation directory."
            time.sleep(2)

    # Build directory
    d = os.path.dirname(os.path.realpath(__file__)) + "/build"
    if os.path.exists(d): 
        dirs.append(d)
    else:
        print "WARNING:\tCannot find build directory in %s. Skipping removal of build directory." % d

    # Data files
    d = os.path.dirname(os.path.realpath(__file__)) + "/src/data"
    if os.path.exists(d):
        files = [d+"/"+f for f in os.listdir(d) if f.endswith(".dat")]
        if len(files) == 0:
            print "WARNING:\tCannot find data files in %s. Skipping removal of data files." % d            
        else:
            for f in files:
                dirs.append(f)
    else:
        print "WARNING:\tCannot find data files in %s. Skipping removal of data directory." % d            

    print ""
        
    if len(dirs) == 0:
        print "WARNING:\tNo files found to delete."
    else:
        while True:
            print('Found the following directories/files:')
            print(str(dirs))
            s = raw_input('Do you really want to delete them?\n[yes/no] ')
            print ""
            if s == "yes":
                # Delete installation files
                for d in dirs:
                    print "Removing %s ..." % d
                    if os.path.isfile(d):
                        os.remove(d)
                    elif os.path.isdir(d):
                        shutil.rmtree("%s" % d, ignore_errors=True)
                    else:
                        print "ERROR: Neither file nor dir: %s" % d
                    print "--> Done"
                break
            elif s == "no":
                print("--> Do not delete files and exit")
                break
            else:
                print("--> Did not understand input \"%s\", please answer either \"yes\" or \"no\"" % s)

