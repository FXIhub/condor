from distutils.core import setup, Extension

nfft_c = Extension('nfft_c', sources = ['nfftmodule.c'], libraries = ['nfft3'], library_dirs = ['/usr/local/lib'],
                   
                   #include_dirs = ["/Users/hantke/.virtualenvs/maxDefault/lib/python2.7/site-packages/numpy/core/include"])
                   include_dirs=['/usr/local/epd-7.2-2-rh5-x86_64/lib'])
                   
                   #include_dirs = ["/Library/Frameworks/Python.framework/Versions/7.2/lib/python2.7/site-packages/numpy/core/include"]
                   #include_dirs = ["/usr/local/epd-7.2-2-rh5-x86_64/lib/python2.7/site-packages/numpy/core/include/"])

setup (name = 'nfft_c',
       version = '1.0',
       description = 'Modules that do some of the nfft functions.',
       ext_modules = [nfft_c])
