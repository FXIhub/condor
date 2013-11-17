import propagator as p
import pylab,os,numpy,sys
import gentools

pdir = os.path.abspath(os.path.dirname(__file__))
odir = pdir+"/example_out/"
os.system("mkdir %s/" % odir)

if len(sys.argv) < 2:
    sample = "all"
else:
    sample = sys.argv[1]


if sample == "all":
    samples = ["icosahedron","sphere","spheroid","uniform_sphere","uniform_spheroid"]
else:
    samples = []
    for i in range(1,len(sys.argv)):
        samples.append(sys.argv[i])

for s in samples:
    print s
    C = gentools.read_configfile(pdir+"/conf/amo55912.conf")
    
    C["sample"] = {}

    if s in ["icosahedron","sphere","spheroid"]:
        C["sample"]["sample_type"] = "map3d"
        C["sample"]["geometry"] = s
        C["sample"]["material_type"] = "virus"
        if s == "spheroid":
            C["sample"]["diameter_c"] = 450E-09
            C["sample"]["diameter_a"] = 250E-09
        else:
            C["sample"]["diameter"] = 450E-09
        I = p.Input(C)
        for i in range(10):
            I.sample.set_random_orientation()
            O = p.Output(I)
            pylab.imsave('%s/example_%s_%i_intensities_poisson.png' % (odir,s,i) ,numpy.log10(numpy.random.poisson(O.get_intensity_pattern())))
            pylab.imsave('%s/example_%s_%i_intensities.png' % (odir,s,i) ,numpy.log10(O.get_intensity_pattern()),vmin=numpy.log10(1.),vmax=numpy.log10(10000.))
            pylab.imsave('%s/example_%s_%i_real_space.png' % (odir,s,i) ,abs(O.get_real_space_image()))

    elif s == "uniform_sphere":
        C["sample"]["sample_type"] = s
        C["sample"]["diameter"] = 450E-09
        C["sample"]["material_type"] = "virus"
        I = p.Input(C)
        O = p.Output(I)
        i = 0
        pylab.imsave('%s/example_%s_%i_intensities_poisson.png' % (odir,s,i) ,numpy.log10(numpy.random.poisson(O.get_intensity_pattern())))
        pylab.imsave('%s/example_%s_%i_intensities.png' % (odir,s,i) ,numpy.log10(O.get_intensity_pattern()),vmin=numpy.log10(1.),vmax=numpy.log10(10000.))
        pylab.imsave('example_%s_%i_real_space.png' % (odir,s,i) ,abs(O.get_real_space_image()))
        
    elif s == "uniform_spheroid":
        C["sample"]["sample_type"] = s
        C["sample"]["diameter_c"] = 450E-09
        C["sample"]["diameter_a"] = 250E-09
        C["sample"]["theta"] = 1.
        C["sample"]["phi"] = 1.
        C["sample"]["material_type"] = "virus"
        I = p.Input(C)
        O = p.Output(I)
        i = 0
        pylab.imsave('%s/example_%s_%i_intensities_poisson.png' % (odir,s,i) ,numpy.log10(numpy.random.poisson(O.get_intensity_pattern())))
        pylab.imsave('%s/example_%s_%i_intensities.png' % (odir,s,i) ,numpy.log10(O.get_intensity_pattern()),vmin=numpy.log10(1.),vmax=numpy.log10(10000.))
        pylab.imsave('%s/example_%s_%i_real_space.png' % (odir,s,i) ,abs(O.get_real_space_image()))

    else:
        print "ERROR: INVALID SAMPLE: %s" % s
