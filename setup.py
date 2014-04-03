# Installation of propagator
import sys, os, fileinput
import constants_data.fetchsf as sf

# Check if python_tools are installed
# try:
#     print "Checking if python_tools are installed..."
#     import imgtools, gentools, cxitools
#     print "Necessary python tools are installed."
# except:
#     print "ERROR: Cannot import python_tools. Please install Max' python_tools and add them to your PYTHONPATH before you proceed."
#     print "Installation of propagator failed."
#     quit(0)

# Scattering factors from the Henke tables and atomic masses 
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
