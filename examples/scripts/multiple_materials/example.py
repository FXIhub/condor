import numpy
import condor

# Core
atomic_composition_core = {"Xe": 1.}
massdensity_core = 3640.
diameter_core = 50E-9

# Shell
atomic_composition_shell = {"H": 2., "O": 1.}
massdensity_shell = 1000.
thickness_shell = 30E-9

S = condor.Source(wavelength=1E-9, pulse_energy=1E-3, focus_diameter=2E-6)
D = condor.Detector(pixel_size=75E-6, nx=1024, ny=1024, distance=0.15)

# Build 3D map for particle model
dx = 0.25*D.get_resolution_element_x(S.photon.get_wavelength())
N = int(numpy.round((1.2*diameter_core+2*thickness_shell)/dx))
assert (dx*N) > (1.1*diameter_core+2*thickness_shell)
map3d = numpy.zeros(shape=(2,N,N,N))
X,Y,Z = numpy.mgrid[:N,:N,:N]
R = numpy.sqrt((X-(N-1)/2.)**2+(Y-(N-1)/2.)**2+(Z-(N-1)/2.)**2)*dx
(map3d[0])[R <= diameter_core/2.] = 1.
if thickness_shell > dx:
    (map3d[1])[(R > diameter_core/2.)*(R <= (diameter_core+thickness_shell*2)/2.)] = 1.

# Initialise particle class instance
P = condor.ParticleMap(geometry='custom',
                       material_type=['custom','custom'], 
                       dx=dx, map3d=map3d,
                       atomic_composition=[atomic_composition_core, atomic_composition_shell],
                       massdensity=[massdensity_core, massdensity_shell])
    
# Initialise experiment class instance
E = condor.Experiment(source=S, particles={"particle_map": P}, detector=D)

# Calculate diffraction pattern
res = E.propagate()

# Images for plotting
img_intensities = res["entry_1"]["data_1"]["data"]
img_fourier = res["entry_1"]["data_1"]["data_fourier"]
real_space = numpy.fft.fftshift(numpy.fft.fft2(numpy.fft.fftshift(img_fourier)))
