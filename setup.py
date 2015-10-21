#!/usr/bin/env python
#--------------------------------------------------------------------------
# CONDOR
# URL: xfel.icm.uu.se/condor
# Copyright 2014 Max Hantke, Filipe R.N.C. Maia, Tomas Ekeberg
# Condor is distributed under the terms of the GNU General Public License
#--------------------------------------------------------------------------

# Always prefer setuptools over distutils
from setuptools import setup, find_packages, Extension
from setuptools.command.install import install
# To use a consistent encoding
from codecs import open
# Other stuff
import sys, os, fileinput
from textwrap import dedent
import numpy

here = os.path.dirname(os.path.realpath(__file__))

# Pickle tables from text files
sys.path.append(here + "/src/data")
import pickle_tables as pt
# Atomic scattering factors from the Henke tables
# B.L. Henke, E.M. Gullikson, and J.C. Davis. X-ray interactions: photoabsorption, scattering, transmission, and reflection at E=50-30000 eV, Z=1-92
# Atomic Data and Nuclear Data Tables Vol. 54 (no.2), 181-342 (July 1993).
# http://henke.lbl.gov/optical_constants/asf.html
print 'Generate data file with atomic scattering constants...'
pt.pickle_atomic_scattering_factors("./src/data/sf","./src/data")
print 'Done.'
# Standard atomic weights from the IUPAC tables
# Atomic weights of the elements 2011 (IUPAC Technical Report) Michael E. Wieser et al., Pure and Applied Chemistry. Volume 85, Issue 5, Pages 1047-1078
# ISSN (Online) 1365-3075, ISSN (Print) 0033-4545
# DOI: 10.1351/PAC-REP-13-03-02, April 2013
# Data loaded from: http://www.chem.qmul.ac.uk/iupac/AtWt/ table 2 (2015/07/01)
print 'Generate data file with atomic standard weight constants...'
pt.pickle_atomic_standard_weights_and_numbers("./src/data/sw","./src/data")
print 'Done.'

# Create example configuration files by concatenation
print "Concatenating example configuration files..."
def concatenate_files(infilenames,outfilename,extra_newlines=0):
    lines = []
    for infilename in infilenames:
        with open(infilename, "r") as f:
            lines += f.readlines()
        if extra_newlines > 0:
            lines += ["\n"]*extra_newlines
    with open(outfilename, "w") as f:
        f.writelines(lines)
src = here + "/examples/configfile/source.conf"
det = here + "/examples/configfile/detector.conf"
for m in ["particle_sphere","particle_spheroid","particle_map","particle_molecule"]:
    par = here + ("/examples/configfile/%s.conf" % m)
    infilenames = [src, par, det]
    d = here + ("/examples/configfile/%s" % m)
    if not os.path.exists(d):
        os.mkdir(d)
        # Ensure that user can write to the directory if setup is run by root
        os.system("chmod a+rwx %s" % d)
    outfilename = d + "/condor.conf"
    print "Concatenating " + outfilename + "..."
    concatenate_files(infilenames, outfilename, 2)
print "Done."

def make_extension_modules(mode="disable_threads", nfft_library_dirs=[], nfft_include_dirs=[]):

    ext_icosahedron = Extension(
        "condor.utils.icosahedron",
        sources=["src/utils/icosahedron/icosahedronmodule.c"],
        include_dirs=[numpy.get_include()],
    )

    _nfft_libraries = {
        "disable_threads": ["nfft3"],
        "enable_threads": ["nfft3_threads" ,"fftw3_threads" ,"fftw3"]
    }
    _nfft_macros = {
        "disable_threads" : [],
        "enable_threads" : [("ENABLE_THREADS", None)],
    }    
    ext_nfft = Extension(
        "condor.utils.nfft",
        sources=["src/utils/nfft/nfftmodule.c"],
        library_dirs=nfft_library_dirs,
        libraries=_nfft_libraries[mode],
        include_dirs=[numpy.get_include()] + nfft_include_dirs,
        define_macros=_nfft_macros[mode],
    )

    return [ext_icosahedron, ext_nfft]

class InstallCommand(install):
    user_options = install.user_options + [
        ('nfft-include-dir=', None, 'Specify the include directory of the NFFT library.'),
        ('nfft-library-dir=', None, 'Specify the library directory of the NFFT library.'),
        ('enable-threads=', None, 'Enable using threads (requires nfft installation with threads (https://www-user.tu-chemnitz.de/~potts/paper/openmpNFFT.pdf).'),
    ]
    
    def initialize_options(self):
        self.nfft_include_dir = None
        self.nfft_library_dir = None
        self.enable_threads   = False
        install.initialize_options(self)
        
    def run(self):
        self.distribution.ext_modules = make_extension_modules(
            "enable_threads" if self.enable_threads else "disable_threads",
            nfft_library_dirs = [self.nfft_library_dir] if self.nfft_library_dir is not None else [],
            nfft_include_dirs = [self.nfft_include_dir] if self.nfft_include_dir is not None else [],
        )
        install.run(self)

# Get the long description from the README file
with open(os.path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()
    
setup(
    cmdclass={
        'install': InstallCommand,
    },
    
    name='condor',

    # Versions should comply with PEP440.  For a discussion on single-sourcing
    # the version across setup.py and the project code, see
    # https://packaging.python.org/en/latest/single_source_version.html
    version='1.0.0',

    description='Condor: Simulation of single particle X-ray diffraction patterns',
    long_description=long_description,

    # The project's main homepage.
    url='https://github.com/mhantke/condor',

    # Author details
    author='Hantke, Max Felix',
    author_email='hantke@xray.bmc.uu.se',

    # Choose your license
    license='GPL',

    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 3 - Alpha',

        # Indicate who your project is intended for
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Physics',

        # Pick your license as you wish (should match "license" above)
        'License :: OSI Approved :: GNU General Public License (GPL)',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        #'Programming Language :: Python :: 3',
        #'Programming Language :: Python :: 3.2',
        #'Programming Language :: Python :: 3.3',
        #'Programming Language :: Python :: 3.4',
    ],

    # What does your project relate to?
    keywords='X-ray diffraction single particle',

    # You can just specify the packages manually here if your project is
    # simple. Or you can use find_packages().
    packages = ['condor', 'condor.utils', 'condor.particle', 'condor.scripts'],
    package_dir = {'condor':'src'},
    
    # Alternatively, if you want to distribute just a my_module.py, uncomment 
    # this:
    #   py_modules=["my_module"],

    # List run-time dependencies here.  These will be installed by pip when
    # your project is installed. For an analysis of "install_requires" vs pip's
    # requirements files see:
    # https://packaging.python.org/en/latest/requirements.html
    install_requires=['numpy', 'scipy', 'h5py'],

    # List additional groups of dependencies here (e.g. development
    # dependencies). You can install these using the following syntax,
    # for example:
    # $ pip install -e .[dev,test]
    #extras_require={
    #    'dev': ['check-manifest'],
    #    'test': ['coverage'],
    #},

    # If there are data files included in your packages that need to be
    # installed, specify them here.  If using Python 2.6 or less, then these
    # have to be included in MANIFEST.in as well.
    package_data={
        'condor': ['data/sf.dat',
                   'data/sw.dat',
                   'data/z.dat',
                   'data/condor.conf',
                   'data/DNA.pdb']
    },

    # Although 'package_data' is the preferred approach, in some case you may
    # need to place data files outside of your packages. See:
    # http://docs.python.org/3.4/distutils/setupscript.html#installing-additional-files # noqa
    # In this case, 'data_file' will be installed into '<sys.prefix>/my_data'
    #data_files=[('my_data', ['data/data_file'])],

    # To provide executable scripts, use entry points in preference to the
    # "scripts" keyword. Entry points provide cross-platform support and allow
    # pip to create the appropriate form of executable for the target platform.
    entry_points={
        'console_scripts': [
            'condor=condor.scripts.condor_script:main',
        ],
    },

    #ext_modules=make_extension_modules(),

)
