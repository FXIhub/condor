import numpy
import os,time
this_dir = os.path.dirname(os.path.realpath(__file__))

try:
    import matplotlib.pyplot as pypl
    plotting = True
except:
    plotting = False

import condor

import logging
logger = logging.getLogger("condor")
#logger.setLevel("DEBUG")
logger.setLevel("WARNING")
#logger.setLevel("INFO")

if False:
    N = 1
    extrinsic_quaternion = numpy.array([ 0.05419534,  0.29519801, -0.91633864,  0.26503677])
    rotation_formalism= "quaternion"
    rotation_values = [extrinsic_quaternion]
else:
    N = 2
    rotation_formalism="random"
    rotation_values = None

# Source
src = condor.Source(wavelength=1.0E-9, pulse_energy=1E-3, focus_diameter=1E-6)
# Detector
det = condor.Detector(distance=1.0, pixel_size=300E-6, nx=256, ny=256)
# Map
#print "Simulating map"
par = condor.ParticleMap(diameter=600E-9, material_type="cell", geometry="custom",
                         map3d_filename="../../map3d.h5", map3d_dataset="data", dx=5E-9,
                         rotation_formalism=rotation_formalism, rotation_values=rotation_values)
s = "particle_map"
E = condor.Experiment(src, {s : par}, det)

W = condor.utils.cxiwriter.CXIWriter("./condor.cxi")
for i in range(N):
    t = time.time()
    res = E.propagate()
    #print time.time()-t
    if plotting:
        real_space = numpy.fft.fftshift(numpy.fft.ifftn(res["entry_1"]["data_1"]["data_fourier"]))
        pypl.imsave(this_dir + "/simple_test_%s_%i.png" % (s,i), numpy.log10(res["entry_1"]["data_1"]["data"]))
        pypl.imsave(this_dir + "/simple_test_%s_%i_rs.png" % (s,i), abs(real_space))
    W.write(res)
W.close()
