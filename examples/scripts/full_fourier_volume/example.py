import os
import condor
this_dir = os.path.dirname(os.path.realpath(__file__))

do_plot = False

if do_plot:
    import matplotlib
    matplotlib.use('TkAgg')
    from matplotlib import pyplot
    from matplotlib.colors import LogNorm

src = condor.Source(wavelength=0.1E-10, pulse_energy=1E-3, focus_diameter=1E-6)
ds = 16
det = condor.Detector(distance=2., pixel_size=110E-6*ds, nx=1024/ds+1, ny=1024/ds+1, solid_angle_correction=False)

psphere = {"particle_sphere": condor.ParticleSphere(diameter=1E-9, material_type="water")}
pmap = {"particle_map": condor.ParticleMap(diameter=1E-9, material_type="water", geometry="icosahedron")}
patoms = {"particle_atoms": condor.ParticleAtoms(pdb_filename="%s/../../DNA.pdb" % this_dir)}
          
particles = [psphere, pmap, patoms]

if do_plot:
    fig, (axs1, axs2) = pyplot.subplots(2, len(particles), figsize=(3*len(particles), 3*2))

for i,par in enumerate(particles): 

    E = condor.Experiment(src, par, det)
    res = E.propagate()
    data = res["entry_1"]["data_1"]["data"]
    if do_plot:
        axs1[i].set_title("2D: " + par.keys()[0])
        lims = (data.min(), data.max())
        axs1[i].imshow(data, norm=LogNorm(lims[0], lims[1]), cmap="gnuplot")
        
    res = E.propagate3d()
    data = res["entry_1"]["data_1"]["data"][1024/ds/2,:,:] 
    if do_plot:
        axs2[i].set_title("3D slice: " + par.keys()[0])
        axs2[i].imshow(data, norm=LogNorm(lims[0], lims[1]), cmap="gnuplot")

if do_plot:
    fig.savefig("2Dvs3D.png", dpi=300)
    pyplot.show()
