#!/usr/bin/env python

#--------------------------------------------------------------------------
# CONDOR
# URL: xfel.icm.uu.se/condor
# Copyright 2014 Max Hantke, Filipe R.N.C. Maia, Tomas Ekeberg
# Condor is distributed under the terms of the GNU General Public License
#--------------------------------------------------------------------------

# Installation of Condor
import sys, os, fileinput
this_dir = os.path.dirname(os.path.realpath(__file__))
print this_dir
sys.path.append(this_dir + "/src/data")
import pickle_tables as pt
from distutils.core import setup, Extension
import numpy.distutils.misc_util

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
src = this_dir + "/examples/configfile/source.conf"
sam = this_dir + "/examples/configfile/sample.conf"
det = this_dir + "/examples/configfile/detector.conf"
for m in ["particle_sphere","particle_spheroid","particle_map","particle_molecule"]:
    par = this_dir + ("/examples/configfile/%s.conf" % m)
    infilenames = [src, sam, par, det]
    d = this_dir + ("/examples/configfile/%s" % m)
    if not os.path.exists(d):
        os.mkdir(d)
        # Ensure that user can write to the directory if setup is run by root
        os.system("chmod a+rwx %s" % d)
    outfilename = d + "/condor.conf"
    print "Concatenating " + outfilename + "..."
    concatenate_files(infilenames, outfilename, 2)
print "Done."

setup(name='condor',
      description='Simulator for diffractive single-particle imaging experiments with X-ray lasers',
      version='1.0',
      author='Max Felix Hantke, Filipe R.N.C. Maia, Tomas Ekeberg',
      author_email='hantke@xray.bmc.uu.se',
      url='http://xfel.icm.uu.se/condor/',
      package_dir={"condor":"src"},
      packages=['condor', 'condor.utils', 'condor.particle'],
      package_data={'condor':['data/sf.dat','data/sw.dat','data/z.dat','data/condor.conf','data/DNA.pdb']},
      scripts=['bin/condor','bin/test_condor'],
      ext_modules=[Extension("condor.utils.icosahedron", sources=["src/utils/icosahedron/icosahedronmodule.c"]),#, extra_compile_args=["-w"]),
                   Extension("condor.utils.nfft", sources=["src/utils/nfft/nfftmodule.c"], libraries=["nfft3"])],#, extra_compile_args=["-w"])],
      include_dirs=numpy.distutils.misc_util.get_numpy_include_dirs()
     )
