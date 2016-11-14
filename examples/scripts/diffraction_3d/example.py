import numpy as np
import scipy.interpolate
import condor

# Number of frames
N = 100

# Dimensions in diffraction space
nx,ny,nz = (100,100,100)

S = condor.Source(wavelength=1E-9, focus_diameter=1E-6, pulse_energy=1E-3)
P = condor.ParticleMap(geometry="icosahedron", diameter=100E-9, material_type="cell", rotation_formalism="random")
D = condor.Detector(pixel_size=1000E-6, distance=0.5, nx=nx, ny=ny)

E = condor.Experiment(source=S, particles={"particle_map": P}, detector=D)

points = []
values = []
for i in range(N):
    res = E.propagate()
    img = res["entry_1"]["data_1"]["data_fourier"]
    qmap = E.get_qmap_from_cache()
    c = 2*np.pi * D.pixel_size / (S.photon.get_wavelength() * D.distance)
    points.append(qmap.reshape((qmap.shape[0]*qmap.shape[1], 3)) / c)
    values.append(img.flatten())
points = np.array(points)
points = points.reshape((points.shape[0]*points.shape[1], 3))
values = np.array(values).flatten()

grid_x, grid_y, grid_z = np.mgrid[0:(nx-1):nx*1j, 0:(ny-1):ny*1j, 0:(nz-1):nz*1j]
grid_x = (grid_x-(nx-1)/2.) 
grid_y = (grid_y-(ny-1)/2.)
grid_z = (grid_z-(nz-1)/2.)
grid = (grid_x, grid_y, grid_z)

# Complex valued 3D diffraction space
img_3d = scipy.interpolate.griddata(points, values, grid, method='linear')

intensities_3d = abs(img_3d)**2
phases_3d = np.angle(img_3d)

# Real space object
tmp = img_3d.copy()
tmp[np.isnan(tmp)] = 0.
tmp = np.fft.fftshift(tmp)
rs_3d = np.fft.fftshift(np.fft.ifftn(tmp))
