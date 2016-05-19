import numpy

try:
    import matplotlib.pyplot as pypl
    plotting = True
except:
    plotting = False
    
import os,time
this_dir = os.path.dirname(os.path.realpath(__file__))

import condor

import logging
logger = logging.getLogger("condor")
#logger.setLevel("DEBUG")
logger.setLevel("WARNING")
#logger.setLevel("INFO")

N = 1
rotation_formalism="random"
rotation_values = None

# Source
src = condor.Source(wavelength=1E-10, pulse_energy=1E-3, focus_diameter=1001E-9)
# Detector
det = condor.Detector(distance=0.2, pixel_size=800E-6, nx=250, ny=250)
# Map
#print "Simulating map"
par = condor.ParticleAtoms(pdb_id="1AKI",
                           rotation_formalism=rotation_formalism, rotation_values=rotation_values)
s = "particle_atoms"
E = condor.Experiment(src, {s : par}, det)

W = condor.utils.cxiwriter.CXIWriter("./condor.cxi")
for i in range(N):
    t = time.time()
    res = E.propagate()
    #print time.time()-t
    if plotting:
        real_space = numpy.fft.fftshift(numpy.fft.ifftn(res["entry_1"]["data_1"]["data_fourier"]))
        pypl.imsave(this_dir + "/%i.png" % (i), numpy.log10(res["entry_1"]["data_1"]["data"]))
        pypl.imsave(this_dir + "/%i_rs.png" % (i), abs(real_space))
    W.write(res)
W.close()

