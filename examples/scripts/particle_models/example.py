import numpy
import matplotlib.pyplot as pypl
import os
this_dir = os.path.dirname(os.path.realpath(__file__))

import condor

import logging
logger = logging.getLogger("condor")
#logger.setLevel("DEBUG")
logger.setLevel("WARNING")
#logger.setLevel("INFO")

# Source
src = condor.Source(wavelength=0.1E-9, pulse_energy=1E-3, focus_diameter=1E-6)
# Sample
sam = condor.Sample()
# Detector
det = condor.Detector(distance=0.05, pixel_size=110E-6, nx=1000, ny=1000)

# Sphere
print "Simulating sphere"
par = condor.ParticleSphere(diameter=1E-9, material_type="water")
s = "sphere"
sam.append_particle(par, s)
E = condor.Experiment(src, sam, det)
res = E.propagate()
real_space = numpy.fft.fftshift(numpy.fft.ifftn(res["entry_1"]["data_1"]["data_fourier"]))
pypl.imsave(this_dir + "/simple_test_%s.png" % s, numpy.log10(res["entry_1"]["data_1"]["data"]))
pypl.imsave(this_dir + "/simple_test_%s_rs.png" % s, abs(real_space))
sam.remove_all_particles()

# Spheroid
print "Simulating spheroid"
par = condor.ParticleSpheroid(diameter=1E-9, material_type="water", flattening=0.6, rotation_formalism="random_z")
s = "spheroid"
sam.append_particle(par, s)
E = condor.Experiment(src, sam, det)
res = E.propagate()
real_space = numpy.fft.fftshift(numpy.fft.ifftn(res["entry_1"]["data_1"]["data_fourier"]))
pypl.imsave(this_dir + "/simple_test_%s.png" % s, numpy.log10(res["entry_1"]["data_1"]["data"]))
pypl.imsave(this_dir + "/simple_test_%s_rs.png" % s, abs(real_space))
sam.remove_all_particles()

# Icosahedron
print "Simulating map"
par = condor.ParticleMap(diameter=1E-9, material_type="water", geometry="icosahedron")
s = "icosahedron"
sam.append_particle(par, s)
E = condor.Experiment(src, sam, det)
res = E.propagate()
real_space = numpy.fft.fftshift(numpy.fft.ifftn(res["entry_1"]["data_1"]["data_fourier"]))
pypl.imsave(this_dir + "/simple_test_%s.png" % s, numpy.log10(res["entry_1"]["data_1"]["data"]))
pypl.imsave(this_dir + "/simple_test_%s_rs.png" % s, abs(real_space))
sam.remove_all_particles()

# Molecule
print "Simulating molecule"
par = condor.ParticleMolecule(pdb_filename="../../DNA.pdb")
s = "molecule"
sam.append_particle(par, s)
E = condor.Experiment(src, sam, det)
res = E.propagate()
real_space = numpy.fft.fftshift(numpy.fft.ifftn(res["entry_1"]["data_1"]["data_fourier"]))
pypl.imsave(this_dir + "/simple_test_%s.png" % s, numpy.log10(res["entry_1"]["data_1"]["data"]))
pypl.imsave(this_dir + "/simple_test_%s_rs.png" % s, abs(real_space))
sam.remove_all_particles()
