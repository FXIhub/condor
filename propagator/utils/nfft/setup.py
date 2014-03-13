from distutils.core import setup, Extension
import numpy.distutils.misc_util
import socket,os
domainname = socket.getfqdn().split(".")[-1]
print domainname
if domainname == "davinci":
    include_dirs = ["/davinci/lib"]
    library_dirs = ["/davinci/lib"]
elif domainname == "gauguin":
    include_dirs = ['/usr/local/epd-7.2-2-rh5-x86_64/lib']
    library_dirs = ['/usr/local/lib']
else:
    include_dirs = numpy.distutils.misc_util.get_numpy_include_dirs()
    #include_dirs = ["/Users/hantke/.virtualenvs/maxDefault/lib/python2.7/site-packages/numpy/core/include","/Library/Frameworks/Python.framework/Versions/7.2/lib/python2.7/site-packages/numpy/core/include","/usr/local/epd-7.2-2-rh5-x86_64/lib"]
    library_dirs = ['/usr/local/lib']
    
if "nfft.so" in os.listdir("."): os.system("rm nfft.so")
if "build" in os.listdir("."): os.system("rm -r build")

nfft = Extension('nfft',
                 sources = ['nfftmodule.c'],
                 libraries = ['nfft3'],
                 library_dirs = library_dirs,
                 include_dirs = include_dirs,
                 extra_compile_args=['-std=c99'])

setup (name = 'nfft',
       version = '1.0',
       description = 'Modules that do some of the nfft functions.',
       ext_modules = [nfft])

sofile = "build/"
sofile += filter(lambda x: x.find("lib.") != -1,os.listdir("build"))[0]+"/"
sofile += "nfft.so"
os.system("ln -s %s" % sofile)
