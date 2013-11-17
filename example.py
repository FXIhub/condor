import propagator as p
import pylab,os,numpy,sys
import gentools

if len(sys.argv) == 1:
    sample = "icosahedron"
else:
    sample = sys.argv[1]

print sample

pdir = os.path.abspath(os.path.dirname(__file__))
C = gentools.read_configfile(pdir+"/conf/amo55912.conf")

C["sample"] = {}
if sample == "icosahedron":
    C["sample"]["sample_type"] = "map3d"
    C["sample"]["geometry"] = "sphere"
    C["sample"]["material_type"] = "virus"
    C["sample"]["diameter"] = 450E-09
    I = p.Input(C)
    for i in range(1):
        #I.sample.set_random_orientation()
        O = p.Output(I)
        print O.get_intensity_pattern()
        #pylab.imsave('example_%i_intensities_poisson.png' % i ,numpy.log10(numpy.random.poisson(O.get_intensity_pattern())))
        pylab.imsave('example_%i_intensities.png' % i ,numpy.log10(O.get_intensity_pattern()),vmin=numpy.log10(1.),vmax=numpy.log10(10000.))
        pylab.imsave('example_%i_real_space.png' % i ,abs(O.get_real_space_image()))
        print "%e" % O.get_intensity_pattern().sum()

elif sample == "sphere":
    C["sample"]["sample_type"] = "uniform_sphere"
    C["sample"]["diameter"] = 450E-09
    C["sample"]["material_type"] = "virus"
    I = p.Input(C)
    O = p.Output(I)
    print "%e" % O.get_intensity_pattern().sum()
    pylab.imsave('example_sphere_intensities_poisson.png' ,numpy.log10(numpy.random.poisson(O.get_intensity_pattern())))
    pylab.imsave('example_sphere_intensities.png' ,numpy.log10(O.get_intensity_pattern()),vmin=numpy.log10(1.),vmax=numpy.log10(10000.))
    pylab.imsave('example_sphere_real_space.png' ,abs(O.get_real_space_image()))
