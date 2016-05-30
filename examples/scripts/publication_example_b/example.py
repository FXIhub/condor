import numpy
import h5py

try:
    import matplotlib.pyplot as pypl
    plotting = True
except:
    plotting = False

import os
this_dir = os.path.dirname(os.path.realpath(__file__))

import condor

import logging
logger = logging.getLogger("condor")
#logger.setLevel("DEBUG")
logger.setLevel("WARNING")
#logger.setLevel("INFO")

wavelength = 1E-10

# Get map with values between 0 and 1
mat = condor.utils.material.AtomDensityMaterial(material_type="custom", massdensity=1340., atomic_composition={"C":332652, "H":492388, "N":98245, "O":131196, "P":7501, "S":2340})
raw_map3d, dx = condor.utils.emdio.fetch_map(1144)
ed_water = condor.utils.material.AtomDensityMaterial(material_type="water").get_electron_density()
ed_particle = mat.get_electron_density()
bin_map3d = condor.utils.emdio.preproc_map_auto(raw_map3d, ed_water=ed_water, ed_particle=ed_particle)
# Scale map to refractive index
dn = mat.get_dn(photon_wavelength=wavelength)
dn_map3d = bin_map3d * dn

# Save map
with h5py.File("emd_1144.h5", "w") as f:
    f["/dn_map"] = dn_map3d
    f["/dx"] = dx
    f["/wavelength"] = wavelength

# Source
src = condor.Source(wavelength=wavelength, pulse_energy=1E-3, focus_diameter=0.1E-6)
# Detector
det = condor.Detector(distance=0.4, pixel_size=125E-6, nx=512, ny=512)
# Particle
par = condor.ParticleMap(diameter=1E-9, geometry="custom", map3d=dn_map3d, dx=dx, rotation_formalism="quaternion", rotation_values=[0.01, 0.69, 0.69, -0.22])

E = condor.Experiment(src, {"particle_map_1144" : par}, det)
res = E.propagate()

W = condor.utils.cxiwriter.CXIWriter("./condor.cxi")
W.write(res)
W.close()

if plotting:
    real_space = numpy.fft.fftshift(numpy.fft.ifftn(res["entry_1"]["data_1"]["data_fourier"]))
    pypl.imsave(this_dir + "/intensities.png", numpy.log10(res["entry_1"]["data_1"]["data"]))
    pypl.imsave(this_dir + "/real_space.png", abs(real_space), cmap="binary")
