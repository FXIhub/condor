import sys,numpy,pylab,multiprocessing
import config
sys.path.append(config.PROPAGATOR_DIR+"/utils/nfft")
import nfft
import time
import logging
logger = logging.getLogger("Propagator")

#def generate_phase_ramp(coordinates,phi,theta,psi):
#    M = pylab.array([[pylab.cos(theta)*pylab.cos(psi),
#                      -pylab.cos(phi)*pylab.sin(psi)+pylab.sin(phi)*pylab.sin(theta)*pylab.cos(psi),
#                      pylab.sin(phi)*pylab.sin(psi)+pylab.cos(phi)*pylab.sin(theta)*pylab.cos(psi)],
#                     [pylab.cos(theta)*pylab.sin(psi),
#                      pylab.cos(phi)*pylab.cos(psi)+pylab.sin(phi)*pylab.sin(theta)*pylab.sin(psi),
#                      -pylab.sin(phi)*pylab.cos(psi)+pylab.cos(phi)*pylab.sin(theta)*pylab.sin(psi)],
#                     [-pylab.sin(theta),
#                       pylab.sin(phi)*pylab.cos(theta),
#                       pylab.cos(phi)*pylab.cos(theta)]])
#    phase_ramp = pylab.zeros(len(coordinates)/3)
#    for i in pylab.arange(0,len(coordinates),3):
#        phase_ramp[i/3] = (coordinates[i]+coordinates[i+1]+coordinates[i+2])*pylab.pi
#    return phase_ramp
    
def arrange_values(values,arrayshape):
    return numpy.reshape(values,arrayshape)

def nfftSingleCore(sample,q):
    invalid_mask = (abs(q)>0.5)
    if (invalid_mask).sum() > 0:
        q[invalid_mask] = 0.
    logger.info("%s invalid pixel positions." % invalid_mask.sum())
    qflat = q.flatten()
    fourierpattern = nfft.nfft3d(qflat,sample)
    fourierpattern = arrange_values(fourierpattern,(q.shape[0],q.shape[1]))
    if (invalid_mask).sum() > 0:
        fourierpattern[numpy.any(invalid_mask,2)] = numpy.nan
    return fourierpattern

def tester():
    s = generate_sample("cube",100)
    result = nfftSingleCore({"y_min":-0.25,
                             "y_steps":500,
                             "x_min":-0.25,
                             "x_steps":200,
                             "stepsize" : 0.01,
                             "phi" : 0.0,
                             "theta" : 0.0,
                             "psi" : 0.0},
                            s)
    print result
    pylab.imsave("testpattern.png",result)

def nfftXCore(object3d,q,N_processes=None):
    #object3d = numpy.complex128(object3d)
    if abs(q).max() > 0.5:
        print "ERROR: nfft accepts only a frequency range from -0.5 to 0.5."
        return
    if N_processes == None:
        N_processes = multiprocessing.cpu_count()
    q_fragments = []
    for n in numpy.arange(0,N_processes,1):
        if n < (N_processes-1): 
            q_fragments.append(q[round(n*q.shape[0]/(1.0*N_processes)):round((n+1)*q.shape[0]/(1.0*N_processes)),:])
        else:
            q_fragments.append(q[round(n*q.shape[0]/(1.0*N_processes)):,:])
    pool = multiprocessing.Pool()
    results = []
    for n in range(0,N_processes):
        config.OUT.write("Start process %i\n" % n)
        results.append(pool.apply_async(nfftSingleCore, (object3d,q_fragments[n])))
    pool.close()
    pool.join()
    A = numpy.zeros(q.shape[0]*q.shape[1],dtype="complex128")
    offset = 0
    for n in range(0,N_processes):
        config.OUT.write("Collect results from process %i\n" % n)
        r = (results[n].get()).flatten()
        A[offset:offset+len(r)] = r[:]
        offset+=len(r)
    return A.reshape((q.shape[0],q.shape[1]))

# For debugging purposes
def generate_cube(size):
    return numpy.ones(shape=(size,size,size))     

def generate_sphere(size):
    x,y,z = numpy.mgrid[0:size,0:size,0:size]
    x -= (size-1)/2.0 
    y -= (size-1)/2.0
    z -= (size-1)/2.0
    r = numpy.sqrt(x**2+y**2+z**2)
    s = numpy.ones(shape=(size,size,size))     
    s[r>size/2.0] = 0.0
    return s

def generate_sample(geometry,size):
    if geometry == "cube":
        return generate_cube(size)
    elif geometry == "sphere":
        return generate_sphere(size)
        
