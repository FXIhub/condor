from distutils.core import setup, Extension

nfft_c = Extension('nfft_c', sources = ['nfftmodule.c'], libraries = ['nfft3'], library_dirs = ['/usr/local/lib'],
                   include_dirs = ["/usr/lib/python2.6/dist-packages/numpy/core/include"])

setup (name = 'nfft_c',
       version = '1.0',
       description = 'Modules that do some of the nfft functions.',
       ext_modules = [nfft_c])
