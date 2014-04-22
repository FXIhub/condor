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
