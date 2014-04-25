# ----------------------------------------------------------------------------------------------------- 
# PENGUIN 
# Simulator for diffractive single-particle imaging experiments with X-ray lasers
# http://xfel.icm.uu.se/penguin/
# ----------------------------------------------------------------------------------------------------- 
# Copyright 2014 Max Hantke, Filipe R.N.C. Maia, Tomas Ekeberg
# Penguin is distributed under the terms of the GNU General Public License
# ----------------------------------------------------------------------------------------------------- 
# General note:
#  All variables are in SI units by default. Exceptions explicit by variable name.
# ----------------------------------------------------------------------------------------------------- 

import os
from distutils.core import setup, Extension
import numpy.distutils.misc_util

if "nfft.so" in os.listdir("."): os.system("rm nfft.so")
if "build" in os.listdir("."): os.system("rm -r build")

ext = Extension("icosahedron", sources=["icosahedronmodule.c"],
                extra_compile_args=['-std=c99'])
setup(ext_modules=[ext], include_dirs=numpy.distutils.misc_util.get_numpy_include_dirs())

sofile = "build/"
sofile += filter(lambda x: x.find("lib.") != -1,os.listdir("build"))[0]+"/"
sofile += "icosahedron.so"
os.system("ln -s %s" % sofile)
