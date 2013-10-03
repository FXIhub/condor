# Installation of propagator
import sys, os, fileinput
import constants_data.fetchsf as sf

# Scattering factors and massesfile from the Henke tables
print 'Loading scattering constants...'
sf.generate_datafile("constants_data/sf",".")
print 'Done.'

print 'Wrapping NFFT...'
pdir = os.path.dirname(os.path.realpath(__file__))
os.chdir("%s/utils/nfft" % pdir)
os.system("python setup.py build")
os.chdir(pdir)
print 'Done.'

in_path = False
if "PYTHONPATH" in os.environ.keys():
    if os.environ["PYTHONPATH"].find(pdir) != -1:
        in_path = True
if not in_path:
    print "For your convenience add \":%s\" to your environmental variable PYTHONPATH." % pdir

print "Installation of propagator ended successfully!" 
