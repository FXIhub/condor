from distutils.core import setup, Extension
import numpy.distutils.misc_util 

# ext = Extension("nfft", sources=["nfftmodule.c"], libraries=["nfft3"], library_dirs=["/usr/local/lib"], include_dirs=["/usr/local/include", numpy.distutils.misc_util.get_numpy_include_dirs()])

# setup(ext_modules=[ext])

#ext = Extension("nfft", sources=["nfftmodule.c", "nfftclassmodule.c"], libraries=["nfft3"])
ext = Extension("nfft", sources=["nfftmodule.c"], libraries=["nfft3"])
setup(ext_modules=[ext], include_dirs=numpy.distutils.misc_util.get_numpy_include_dirs())
