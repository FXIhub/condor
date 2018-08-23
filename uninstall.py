#!/usr/bin/env python

# -----------------------------------------------------------------------------------------------------
# CONDOR
# Simulator for diffractive single-particle imaging experiments with X-ray lasers
# http://xfel.icm.uu.se/condor/
# -----------------------------------------------------------------------------------------------------
# Copyright 2016 Max Hantke, Filipe R.N.C. Maia, Tomas Ekeberg
# Condor is distributed under the terms of the BSD 2-Clause License
# -----------------------------------------------------------------------------------------------------
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# 
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# -----------------------------------------------------------------------------------------------------
# General note:
# All variables are in SI units by default. Exceptions explicit by variable name.
# -----------------------------------------------------------------------------------------------------

from __future__ import print_function, absolute_import # Compatibility with python 2 and 3
import os,sys
import time
import shutil

import argparse, subprocess

from builtins import input
from distutils.core import setup, Extension
from distutils.spawn import find_executable
from distutils.sysconfig import get_python_lib

parser = argparse.ArgumentParser(description='Deinstallation script for Condor')
parser.add_argument('-p', '--packages-dir', metavar='packages', type=str, 
                    help="packages directory where condor was installed")
args = parser.parse_args()


if __name__ == "__main__":
   
    # Check existence of installation directory
    if args.packages_dir:
        plibs = [get_python_lib(prefix=args.prefix)]
    else:
        plibs = [get_python_lib()]
    dirs = []
    for plib in plibs:
        if os.path.exists(plib):
            dirs += [(plib + "/" + d) for d in os.listdir(plib) if d.startswith("condor")]
        
    if len(dirs) == 0:
        try:
            import condor
            dirs = [os.path.dirname(os.path.realpath(condor.__file__))]
        except:
            print("WARNING:\tCannot find or cannot access condor installation directory")
            print("\t\tUse the --packages-dir argument if your installation is not in a standard location.")

    # Build directory
    d = os.path.dirname(os.path.realpath(__file__)) + "/build"
    if os.path.exists(d): 
        dirs.append(d)
    else:
        print("WARNING:\tCannot find build directory in %s. Skipping removal of build directory." % d)

    # Data files
    d = os.path.dirname(os.path.realpath(__file__)) + "/src/data"
    if os.path.exists(d):
        files = [d+"/"+f for f in os.listdir(d) if f.endswith(".dat")]
        if len(files) == 0:
            print("WARNING:\tCannot find data files in %s. Skipping removal of data files." % d)
        else:
            for f in files:
                dirs.append(f)
    else:
        print("WARNING:\tCannot find data files in %s. Skipping removal of data directory." % d)
    print("")
        
    if len(dirs) == 0:
        print("WARNING:\tNo files found to delete.")
    else:
        while True:
            print('Found the following directories/files:')
            print(str(dirs))
            s = input('Do you really want to delete them?\n[yes/no] ')
            print("")
            if s == "yes":
                # Delete installation files
                for d in dirs:
                    print("Removing %s ..." % d)
                    if os.path.isfile(d):
                        os.remove(d)
                    elif os.path.isdir(d):
                        shutil.rmtree("%s" % d, ignore_errors=True)
                    else:
                        print("ERROR: Neither file nor dir: %s" % d)
                    print("--> Done")
                break
            elif s == "no":
                print("--> Do not delete files and exit")
                break
            else:
                print("--> Did not understand input \"%s\", please answer either \"yes\" or \"no\"" % s)

