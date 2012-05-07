import sys,pylab,multiprocessing
import config
sys.path.append(config.PROPAGATOR_DIR+"/utils/nfft")
import nfft
import time

#def generate_fourier_coordinates(y_min,y_steps,x_min,x_steps,stepsize):
#    X,Y = pylab.meshgrid((pylab.arange(0,x_steps,1)*stepsize+x_min),(pylab.arange(0,y_steps,1)*stepsize+y_min))
#    coords = pylab.zeros(shape=(X.shape[0]*X.shape[1],3))
#    coords[:,1] = Y.flatten()[:]
#    coords[:,2] = X.flatten()[:]
#    return coords.flatten()

#def rotate_coordinates(coordinates,phi,theta,psi):
    # Turn grid and compute phase ramp (correction for translation)
#    M = pylab.array([[pylab.cos(theta)*pylab.cos(psi),
#                      -pylab.cos(phi)*pylab.sin(psi)+pylab.sin(phi)*pylab.sin(theta)*pylab.cos(psi),
#                      pylab.sin(phi)*pylab.sin(psi)+pylab.cos(phi)*pylab.sin(theta)*pylab.cos(psi)],
#                     [pylab.cos(theta)*pylab.sin(psi),
#                      pylab.cos(phi)*pylab.cos(psi)+pylab.sin(phi)*pylab.sin(theta)*pylab.sin(psi),
#                      -pylab.sin(phi)*pylab.cos(psi)+pylab.cos(phi)*pylab.sin(theta)*pylab.sin(psi)],
#                     [-pylab.sin(theta),
#                       pylab.sin(phi)*pylab.cos(theta),
#                       pylab.cos(phi)*pylab.cos(theta)]])
#    for i in pylab.arange(0,len(coordinates),3):
#        [coordinates[i],coordinates[i+1],coordinates[i+2]] = pylab.dot(M,coordinates[i:i+3])
#    return coordinates

def generate_phase_ramp(coordinates,phi,theta,psi):
    M = pylab.array([[pylab.cos(theta)*pylab.cos(psi),
                      -pylab.cos(phi)*pylab.sin(psi)+pylab.sin(phi)*pylab.sin(theta)*pylab.cos(psi),
                      pylab.sin(phi)*pylab.sin(psi)+pylab.cos(phi)*pylab.sin(theta)*pylab.cos(psi)],
                     [pylab.cos(theta)*pylab.sin(psi),
                      pylab.cos(phi)*pylab.cos(psi)+pylab.sin(phi)*pylab.sin(theta)*pylab.sin(psi),
                      -pylab.sin(phi)*pylab.cos(psi)+pylab.cos(phi)*pylab.sin(theta)*pylab.sin(psi)],
                     [-pylab.sin(theta),
                       pylab.sin(phi)*pylab.cos(theta),
                       pylab.cos(phi)*pylab.cos(theta)]])
    phase_ramp = pylab.zeros(len(coordinates)/3)
    for i in pylab.arange(0,len(coordinates),3):
        phase_ramp[i/3] = (coordinates[i]+coordinates[i+1]+coordinates[i+2])*pylab.pi
    return phase_ramp
    
def arrange_values(values,arrayshape):
    return values.reshape(arrayshape)

def nfftSingleCore(sample,q):
    #phase_ramp = generate_phase_ramp(coordinates,params["phi"],params["theta"],params["psi"])
    qflat = q.flatten()
    fourierpattern = nfft.nfft3d(qflat,sample)
    fourierpattern = arrange_values(fourierpattern,(q.shape[0],q.shape[1]))
    #fourierpattern = arrange_values(fourierpattern,(params["y_steps"],params["x_steps"]))
    #phase_ramp = arrange_values(phase_ramp,(params["y_steps"],params["x_steps"]))
    #phases = pylab.log(fourierpattern/abs(fourierpattern)).imag
    #phases -= phase_ramp
    #while len(phases[abs(phases)>pylab.pi]) != 0:
    #    phases[phases<=-pylab.pi] += 2*pylab.pi
    #    phases[phases>pylab.pi] -= 2*pylab.pi
    #return abs(fourierpattern)*pylab.exp(phases*1.0j)
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
    object3d = pylab.complex128(object3d)
    #if len(q.shape) == 3: q = q.reshape((q.shape[0]*q.shape[1],3))
    if abs(q).max() > 0.5:
        print "ERROR: nfft accepts only a frequency range from -0.5 to 0.5."
        return
    if N_processes == None:
        N_processes = multiprocessing.cpu_count()
    q_fragments = []
    for n in pylab.arange(0,N_processes,1):
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
    A = pylab.zeros(q.shape[0]*q.shape[1],dtype="complex128")
    offset = 0
    for n in range(0,N_processes):
        config.OUT.write("Collect results from process %i\n" % n)
        r = (results[n].get()).flatten()
        A[offset:offset+len(r)] = r[:]
        offset+=len(r)
    return A.reshape((q.shape[0],q.shape[1]))

# For debugging purposes
def generate_cube(size):
    return pylab.ones(shape=(size,size,size))     

def generate_sphere(size):
    x,y,z = pylab.mgrid[0:size,0:size,0:size]
    x -= (size-1)/2.0 
    y -= (size-1)/2.0
    z -= (size-1)/2.0
    r = pylab.sqrt(x**2+y**2+z**2)
    s = pylab.ones(shape=(size,size,size))     
    s[r>size/2.0] = 0.0
    return s

def generate_sample(geometry,size):
    if geometry == "cube":
        return generate_cube(size)
    elif geometry == "sphere":
        return generate_sphere(size)
        
