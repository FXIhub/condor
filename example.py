import propagator as p
import pylab

I = p.Input('propagator.conf')
print 'Creating 3d refractive index map of icosahedral virus...'
I.set_sample_virus_map(100E-09)
print 'Done'
print 'Propagation to detector...'
O = p.propagator(I)
print 'Done'
print 'Generation png output'
pylab.imsave('intensities.png',log10(O.get_intensity_pattern()))
print 'Done'
print 'Showing result in new window'
O.plot_pattern(noise='poisson')
print 'Done'
