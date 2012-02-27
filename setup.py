import sys, os, fileinput
import constants_data.fetchsf as sf

print 'Generate file of scattering factors...'
sf.generate_datafile("constants_data/sf",".")
print 'Done.'

print 'Setting up directory...'
PROPAGATOR_DIR = raw_input("Please enter the base path of PROPAGATOR: ")
for line in fileinput.input("constants.py", inplace=1):
    if "PROPAGATOR_DIR =" in line:
        print "PROPAGATOR_DIR = \"%s\"\n" % PROPAGATOR_DIR
    else:
        print line[:-1]
print 'Done.'

print 'Setup ended successfully.'
