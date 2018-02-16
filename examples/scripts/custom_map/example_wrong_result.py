import numpy
import matplotlib.pyplot as pypl
import os
this_dir = os.path.dirname(os.path.realpath(__file__))

import condor

import logging
logger = logging.getLogger("condor")
#logger.setLevel("DEBUG")
#logger.setLevel("WARNING")
logger.setLevel("INFO")

if False:
    N = 1
    extrinsic_quaternion = numpy.array([ 0.05419534,  0.29519801, -0.91633864,  0.26503677])
    rotation_formalism= "quaternion"
    rotation_values = [extrinsic_quaternion]
else:
    N = 2
    rotation_formalism="random"
    rotation_values = None

focus_diameter = 100e-9
intensity = 2*1.02646137e9
pulse_energy = ((focus_diameter**2) * numpy.pi / 4.) * intensity
wavelength=0.2262e-9
pixelsize = 110e-6
distance=2.4
dx = 3.76e-10

# Source
src = condor.Source(wavelength=wavelength, pulse_energy=pulse_energy , focus_diameter=focus_diameter)
# Detector
det = condor.Detector(distance=distance, pixel_size=pixelsize, nx=414, ny=414)
# Map
print("Simulating map")
par = condor.ParticleMap(diameter=40E-9, material_type="poliovirus", geometry="custom",
                         map3d_filename="../../map3d.h5", map3d_dataset="data", dx=dx,
                         rotation_formalism=rotation_formalism, rotation_values=rotation_values)
s = "particle_map"
E = condor.Experiment(src, {s : par}, det)

W = condor.utils.cxiwriter.CXIWriter("./condor.cxi")
for i in range(N):
    res = E.propagate()
    real_space = numpy.fft.fftshift(numpy.fft.ifftn(res["entry_1"]["data_1"]["data_fourier"]))
    pypl.imsave(this_dir + "/simple_test_%s_%i.png" % (s,i), numpy.log10(res["entry_1"]["data_1"]["data"]))
    pypl.imsave(this_dir + "/simple_test_%s_%i_rs.png" % (s,i), abs(real_space))
    W.write(res)
W.close()
