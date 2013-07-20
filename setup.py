import sys, os, fileinput
import constants_data.fetchsf as sf

print 'Generate file of scattering factors...'
sf.generate_datafile("constants_data/sf",".")
print 'Done.'

print 'Set-up nfft wrapper'
pdir = os.path.dirname(os.path.realpath(__file__))
os.chdir("%s/utils/nfft" % pdir)
os.system("python setup.py build")
os.chdir(pdir)
print 'done'
