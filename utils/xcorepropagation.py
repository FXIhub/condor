import sys,pylab,multiprocessing
import constants
sys.path.append(constants.PROPAGATOR_DIR+"/utils/nfft")
import nfft
import time

def generate_fourier_coordinates(y_min,y_steps,x_min,x_steps,stepsize):
    X,Y = pylab.meshgrid((pylab.arange(0,x_steps,1)*stepsize+x_min),(pylab.arange(0,y_steps,1)*stepsize+y_min))
    coords = pylab.zeros(shape=(X.shape[0]*X.shape[1],3))
    coords[:,1] = Y.flatten()[:]
    coords[:,2] = X.flatten()[:]
    return coords.flatten()

def rotate_coordinates(coordinates,phi,theta,psi):
    # Turn grid and compute phase ramp (correction for translation)
    M = pylab.array([[pylab.cos(theta)*pylab.cos(psi),
                      -pylab.cos(phi)*pylab.sin(psi)+pylab.sin(phi)*pylab.sin(theta)*pylab.cos(psi),
                      pylab.sin(phi)*pylab.sin(psi)+pylab.cos(phi)*pylab.sin(theta)*pylab.cos(psi)],
                     [pylab.cos(theta)*pylab.sin(psi),
                      pylab.cos(phi)*pylab.cos(psi)+pylab.sin(phi)*pylab.sin(theta)*pylab.sin(psi),
                      -pylab.sin(phi)*pylab.cos(psi)+pylab.cos(phi)*pylab.sin(theta)*pylab.sin(psi)],
                     [-pylab.sin(theta),
                       pylab.sin(phi)*pylab.cos(theta),
                       pylab.cos(phi)*pylab.cos(theta)]])
    for i in pylab.arange(0,len(coordinates),3):
        [coordinates[i],coordinates[i+1],coordinates[i+2]] = pylab.dot(M,coordinates[i:i+3])
    return coordinates

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

def arrange_values(values,arrayshape):
    return values.reshape(arrayshape)

def nfftSingleCore(sample,params):
    #sample = array(sample_raw,shape=(params["Nz_object"],params["Ny_object"],params["Nx_object"]))
    #print "7 %i" % (psutil.avail_phymem())
    coordinates = generate_fourier_coordinates(params["y_min"],params["y_steps"],params["x_min"],params["x_steps"],params["stepsize"])
    #print "8 %i" % (psutil.avail_phymem())
    #A= coordinates.shape[0]
    coordinates = rotate_coordinates(coordinates,params["phi"],params["theta"],params["psi"])
    #print "9 %i" % (psutil.avail_phymem())
    #B= coordinates.shape[0]
    phase_ramp = generate_phase_ramp(coordinates,params["phi"],params["theta"],params["psi"])
    #print "10 %i" % (psutil.avail_phymem())
    fourierpattern = nfft.nfft3d(coordinates,sample)
    #print "11 %i" % (psutil.avail_phymem())
    #C= values.shape[0]#values.shape[1]
    #D= (params["y_steps"]*params["x_steps"])
    #print "Yeeah"
    #print A/3.0
    #print B/3.0
    #print C
    #print D
    fourierpattern = arrange_values(fourierpattern,(params["y_steps"],params["x_steps"]))
    #print "12 %i" % (psutil.avail_phymem())
    phase_ramp = arrange_values(phase_ramp,(params["y_steps"],params["x_steps"]))
    phases = pylab.log(fourierpattern/abs(fourierpattern)).imag
    phases -= phase_ramp
    while len(phases[abs(phases)>pylab.pi]) != 0:
        phases[phases<=-pylab.pi] += 2*pylab.pi
        phases[phases>pylab.pi] -= 2*pylab.pi
    return abs(fourierpattern)*pylab.exp(phases*1.0j)

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

def nfftXCore(object3d,N,eulerangles=[[0.0,0.0,0.0]],N_processes=None,interval=0.5):
    #object3d = generate_sample("sphere",100)
    object3d = pylab.complex128(object3d)

    if N_processes == None:
        #phymem = psutil.avail_phymem()
        #objectsize = object3d.size*object3d.itemsize
        #print "Mem/Obj %f" % (phymem / (1.0*objectsize))
        N_processes = 2#multiprocessing.cpu_count()

    patterns = []
    stepsize = interval/(1.0*(N-1))
    s_min = -interval/2.0
    s_max = interval/2.0
    ysteps_process = pylab.floor(interval/(N_processes*stepsize))
    [phi,theta,psi] = eulerangles
    processparameters = []
    stepsum = 0

    for n in pylab.arange(0,N_processes,1):
        x_min = s_min
        x_steps = N
        y_min = s_min + n*ysteps_process*stepsize
        if n < (N_processes-1): y_steps = ysteps_process
        else: y_steps = N - ysteps_process*(N_processes-1)
        processparameters.append({"y_min":y_min,
                                  "y_steps":y_steps,
                                  "x_min":x_min,
                                  "x_steps":x_steps,
                                  "stepsize" : stepsize,
                                  "phi" : phi,
                                  "theta" : theta,
                                  "psi" : psi})


    #object3d = object3d.flatten()
    #Object3d = multiprocessing.Array('f',object3d)
    pool = multiprocessing.Pool()
    results = []
    for n in range(0,N_processes):
        #print "4/%i %i" % (n,psutil.avail_phymem())
        results.append(pool.apply_async(nfftSingleCore, (object3d,processparameters[n])))
    #pool.close()
    #print "Available memory: %i bytes" % psutil.avail_phymem()
    while True:
        time.sleep(0.01)
        if psutil.avail_phymem() < 10000000:
            #print "Available memory: %i bytes" % psutil.avail_phymem()
            #print "ERROR: Running out of physical memory. Use less CPUs for propagation."
            pool.terminate()
            #return pylab.zeros((N,N))
        ready_states = []
        for process in results:
            ready_states.append(process.ready())
        #print ready_states
        if pylab.all(ready_states):
            break
    pool.close()
    pool.join()
    A = pylab.zeros(N*N,dtype="complex128")
    offset = 0
    for n in range(0,N_processes):
        r = (results[n].get()).flatten()
        #print r.shape
        #print A.shape
        A[offset:offset+len(r)] = r[:]
        offset+=len(r)

    return A.reshape((N,N))


        
        
