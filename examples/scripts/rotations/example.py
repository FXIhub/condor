import numpy
import matplotlib.pyplot as pypl
import os, shutil
this_dir = os.path.dirname(os.path.realpath(__file__))

import condor

import logging
logger = logging.getLogger('condor')
logger.setLevel("INFO")
#logger.setLevel("DEBUG")

out_dir = this_dir + "/pngs"

if os.path.exists(out_dir):
    shutil.rmtree(out_dir)
os.mkdir(out_dir)

# Source
src = condor.Source(wavelength=0.1E-9, pulse_energy=1E-3, focus_diameter=1E-6)
# Detector
det = condor.Detector(distance=0.5, pixel_size=750E-6, nx=100, ny=100)#, cx=55, cy=55)


#angles_d = numpy.array([0., 22.5, 45.])
angles_d = numpy.array([72.5])

for angle_d in angles_d:

    angle = angle_d/360.*2*numpy.pi
    rotation_axis = numpy.array([1.,1.,0.])/numpy.sqrt(2.)
    quaternion = condor.utils.rotation.quat(angle,rotation_axis[0],rotation_axis[1], rotation_axis[2])
    rotation_values = numpy.array([quaternion])
    rotation_formalism = "quaternion"
    rotation_mode = "extrinsic"
    #rotation_values = None
    #rotation_formalism = "random"
    #rotation_mode = "extrinsic"   
    #rotation_values = None
    #rotation_formalism = None
    #rotation_mode = "extrinsic"

    print "Angle = %.2f degrees" % angle_d

    short_diameter = 25E-9*12/100.
    long_diameter = 2*short_diameter
    spheroid_diameter   = condor.utils.spheroid_diffraction.to_spheroid_diameter(short_diameter/2.,long_diameter/2.)
    spheroid_flattening = condor.utils.spheroid_diffraction.to_spheroid_flattening(short_diameter/2.,long_diameter/2.)
    N_long = 20
    N_short = int(round(short_diameter/long_diameter * N_long))

    # Spheroid
    if True:
        # Ideal spheroid
        print "Simulating spheroid"
        par = condor.ParticleSpheroid(diameter=spheroid_diameter, material_type="water", flattening=spheroid_flattening, rotation_values=rotation_values, rotation_formalism=rotation_formalism, rotation_mode=rotation_mode)
        s = "particle_spheroid"
        E = condor.Experiment(src, {s : par}, det)
        res = E.propagate()
        real_space = numpy.fft.fftshift(numpy.fft.ifftn(res["entry_1"]["data_1"]["data_fourier"]))
        vmin = numpy.log10(res["entry_1"]["data_1"]["data"].max()/10000.)
        pypl.imsave(out_dir + "/%s_%2.2fdeg.png" % (s,angle_d), numpy.log10(res["entry_1"]["data_1"]["data"]), vmin=vmin)
        pypl.imsave(out_dir + "/%s_rs_%2.2fdeg.png" % (s,angle_d), abs(real_space))
    
    if True:
        # Map (spheroid)
        print "Simulating map (spheroid)"
        par = condor.ParticleMap(diameter=spheroid_diameter, material_type="water", flattening=spheroid_flattening, geometry="spheroid", rotation_values=rotation_values, rotation_formalism=rotation_formalism, rotation_mode=rotation_mode)
        s = "particle_map_spheroid"
        E = condor.Experiment(src, {s : par}, det)
        res = E.propagate()
        real_space = numpy.fft.fftshift(numpy.fft.ifftn(res["entry_1"]["data_1"]["data_fourier"]))
        vmin = numpy.log10(res["entry_1"]["data_1"]["data"].max()/10000.)
        pypl.imsave(out_dir + "/%s_%2.2f.png" % (s,angle_d), numpy.log10(res["entry_1"]["data_1"]["data"]), vmin=vmin)
        pypl.imsave(out_dir + "/%s_rs_%2.2f.png" % (s,angle_d), abs(real_space))

    # Box
    if True:
        # Map (box)
        dx = long_diameter/(N_long-1)
        map3d = numpy.zeros(shape=(N_long,N_long,N_long))
        map3d[:N_short,:,:N_short] = 1.
        map3d[N_short:N_short+N_short,:N_short,:N_short] = 1.
        # Map
        print "Simulating map (custom)"
        par = condor.ParticleMap(diameter=long_diameter, material_type="water", geometry="custom", map3d=map3d, dx=dx, rotation_values=rotation_values, rotation_formalism=rotation_formalism, rotation_mode=rotation_mode)
        s = "particle_map_custom"
        E = condor.Experiment(src, {s : par}, det)
        res = E.propagate()
        data_fourier = res["entry_1"]["data_1"]["data_fourier"]
        #data_fourier = abs(data_fourier)*numpy.exp(-1.j*numpy.angle(data_fourier))
        real_space = numpy.fft.fftshift(numpy.fft.ifftn(numpy.fft.fftshift(data_fourier)))
        vmin = numpy.log10(res["entry_1"]["data_1"]["data"].max()/10000.)
        pypl.imsave(out_dir + "/%s_map.png" % (s),map3d.sum(0))
        pypl.imsave(out_dir + "/%s_%2.2f.png" % (s,angle_d), numpy.log10(res["entry_1"]["data_1"]["data"]), vmin=vmin)
        pypl.imsave(out_dir + "/%s_%2.2f_phases.png" % (s,angle_d), numpy.angle(res["entry_1"]["data_1"]["data_fourier"])%(2*numpy.pi))
        pypl.imsave(out_dir + "/%s_rs_%2.2f.png" % (s,angle_d), abs(real_space))

    if True:
        # Atoms (box)
        print "Simulating atoms"
        Z1,Y1,X1 = numpy.meshgrid(numpy.linspace(0, short_diameter, N_short),
                                  numpy.linspace(0, long_diameter,   N_long),
                                  numpy.linspace(0, short_diameter, N_short),
                                  indexing="ij")
        Z2,Y2,X2 = numpy.meshgrid(numpy.linspace(0, short_diameter, N_short) + long_diameter/2.,
                                  numpy.linspace(0, short_diameter, N_short),
                                  numpy.linspace(0, short_diameter, N_short),
                                  indexing="ij")
        Z = numpy.concatenate((Z1.ravel(),Z2.ravel()))
        Y = numpy.concatenate((Y1.ravel(),Y2.ravel()))
        X = numpy.concatenate((X1.ravel(),X2.ravel()))
        proj = numpy.zeros(shape=(N_long,N_long))
        dx = long_diameter/(N_long-1)
        for (x,y,z) in zip(X.ravel(),Y.ravel(),Z.ravel()):
            proj[int(round(y/dx)),int(round(x/dx))] += 1
        pypl.imsave(out_dir + "/%s_proj.png" % (s),proj)
        atomic_positions = numpy.array([[x,y,z] for x,y,z in zip(X.ravel(),Y.ravel(),Z.ravel())])
        atomic_numbers   = numpy.ones(atomic_positions.size/3, dtype=numpy.int16)
        par = condor.ParticleAtoms(atomic_positions=atomic_positions, atomic_numbers=atomic_numbers, rotation_values=rotation_values, rotation_formalism=rotation_formalism, rotation_mode=rotation_mode)
        s = "particle_atoms"
        E = condor.Experiment(src, {s : par}, det)
        res = E.propagate()
        real_space = numpy.fft.fftshift(numpy.fft.ifftn(numpy.fft.fftshift(res["entry_1"]["data_1"]["data_fourier"])))
        fourier_space = res["entry_1"]["data_1"]["data_fourier"]
        vmin = numpy.log10(res["entry_1"]["data_1"]["data"].max()/10000.)
        pypl.imsave(out_dir + "/%s_%2.2f.png" % (s,angle_d), numpy.log10(res["entry_1"]["data_1"]["data"]), vmin=vmin)
        pypl.imsave(out_dir + "/%s_%2.2f_phases.png" % (s,angle_d), numpy.angle(fourier_space)%(2*numpy.pi))
        pypl.imsave(out_dir + "/%s_rs_%2.2f.png" % (s,angle_d), abs(real_space))
    
