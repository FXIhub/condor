from distutils.core import setup, Extension
import numpy.distutils.misc_util 

ext = Extension("nfft", sources=["nfftmodule.c"], libraries=["nfft3"])
setup(name="condor_nfft",ext_modules=[ext], include_dirs=numpy.distutils.misc_util.get_numpy_include_dirs())
