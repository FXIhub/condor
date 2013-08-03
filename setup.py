#!/usr/bin/env python

from distutils.core import setup,Extension
import os

# clean up
if "build" in os.listdir("."): os.system("rm -r build")

print "Installing NFFT wrapper"

# determine numpy include directory
try:
    from numpy.distutils.misc_util import get_numpy_include_dirs
    include_dirs = get_numpy_include_dirs()
    N = len(include_dirs)
    for i in range(N): include_dirs.append(include_dirs[i]+"/numpy")
    print "Found numpy include directories: ",include_dirs
except:
    print "Could not find numpy include directory."
    include_dirs = []
    include_dirs.append(raw_input("Please enter the numply include directory [typically ending with \'site-packages/numpy/core/include/numpy\']: "))

# determine nfft3 library directory
libfound = False
libpaths = os.environ["LD_LIBRARY_PATH"].split(":")
libpaths.append("/lib")
libpaths.append("/usr/lib")
for p in libpaths:
    if "libnfft3.so" in os.listdir(p):
        print "Found libnfft3.so in %s" % p
        libfound = True
        break
if libfound:
    library_dirs = [p]
else:
    print "Could not find nfft3 library."
    library_dirs = []
    library_dirs.append(raw_input("Please enter the nfft3 library directory: "))

nfft_c = Extension('nfft_c', sources = ['src/nfft/nfftmodule.c'], libraries = ['nfft3'], include_dirs = include_dirs)
setup (name = 'nfft_c',
       version = '1.0',
       description = 'Modules that do some of the nfft functions.',
       ext_modules = [nfft_c])

print "Installing main propagator code"

files = ["utils/*"]

print "Generate data-file of scattering factors from Henke tables."
import const.fetchsf as sf
sf.generate_datafile("const/sf",".")

setup(name = "propagator",
      version = "1.0",
      description = "Data simulatior for single particle X-ray diffraction experiments",
      author = "Max Felix Hantke",
      author_email = "maxhantke@gmail.com",
      url = "https://github.com/mhantke/propagator",
      packages = ['propagator'],
      package_dir = {"propagator" : "src/propagator"},
      package_data = {'propagator' : files },
      data_files = [("henke_tables",["elements.dat"])])

print "Installation of propagator finished."
