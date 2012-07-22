import sys, os, fileinput
import constants_data.fetchsf as sf

print 'Generate file of scattering factors...'
sf.generate_datafile("constants_data/sf",".")
print 'Done.'

print 'Setting up directory...'
PROPAGATOR_DIR = raw_input("Please enter the base path of PROPAGATOR: ")
if PROPAGATOR_DIR[-1] == '/': PROPAGATOR_DIR = PROPAGATOR_DIR[:-1]
f = open('_config.py','w')
f.writelines(["# Personal configuration file for propagator\n",
              "PROPAGATOR_DIR = \"%s\"\n" % PROPAGATOR_DIR])
f.close()
print 'Done.'

print 'Setup ended successfully.'
