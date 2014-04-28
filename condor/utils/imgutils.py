# ----------------------------------------------------------------------------------------------------- 
# CONDOR 
# Simulator for diffractive single-particle imaging experiments with X-ray lasers
# http://xfel.icm.uu.se/condor/
# ----------------------------------------------------------------------------------------------------- 
# Copyright 2014 Max Hantke, Filipe R.N.C. Maia, Tomas Ekeberg
# Condor is distributed under the terms of the GNU General Public License
# ----------------------------------------------------------------------------------------------------- 
# General note:
#  All variables are in SI units by default. Exceptions explicit by variable name.
# ----------------------------------------------------------------------------------------------------- 

import os,re,sys,h5py,numpy,time
import condor.config
import condortools

def radial_pixel_average(image,**kargs):
    if 'cx' in kargs: cx = kargs['cx']
    else: cx = (image.shape[1]-1)/2.0
    if 'cy' in kargs: cy = kargs['cy'] 
    else: cy = (image.shape[0]-1)/2.0
    x = numpy.arange(0,image.shape[1],1.0)-cx
    y = numpy.arange(0,image.shape[1],1.0)-cy
    X,Y = numpy.meshgrid(x,y)
    R = numpy.sqrt(X**2+Y**2)
    R = R.round()
    R[image==numpy.Inf] = -1
    radii = numpy.arange(R.min(),R.max()+1,1)
    values = numpy.zeros_like(radii)
    for i in range(0,len(radii)):
        values[i] = image[R==radii[i]].mean()
    if 'rout' in kargs: return numpy.array([radii,values])
    else:return values
        
def radial_pixel_sum(image,**kargs):
    if 'cx' in kargs: cx = kargs['cx']
    else: cx = (image.shape[1]-1)/2.0
    if 'cy' in kargs: cy = kargs['cy'] 
    else: cy = (image.shape[0]-1)/2.0
    x = numpy.arange(0,image.shape[1],1.0)-cx
    y = numpy.arange(0,image.shape[1],1.0)-cy
    X,Y = numpy.meshgrid(x,y)
    R = numpy.sqrt(X**2+Y**2)
    R = R.round()
    radii = numpy.arange(R.min(),R.max()+1,1)
    values = numpy.zeros_like(radii)
    for i in range(0,len(radii)):
        values[i] = image[R==radii[i]].sum()
    if 'rout' in kargs: return numpy.array([radii,values])
    else:return values

def cartesian_to_polar(cartesian_pattern,N_theta,x_center=None,y_center=None):
    Nx = cartesian_pattern.shape[1]
    Ny = cartesian_pattern.shape[0]
    R = int(min([Nx,Ny])/2.0-1)
    polar_pattern = numpy.zeros(shape=(R,N_theta))
    if not x_center:
        x_center = Nx/2.0-0.5
    if not y_center:
        y_center = Ny/2.0-0.5
    for i_theta in range(0,N_theta):
        for r in range(0,R):
            theta = 2*numpy.pi*i_theta/float(N_theta)
            x = x_center + r * numpy.sin(theta)
            y = y_center + r * numpy.cos(theta)
            # bilinear interpolation
            x1 = int(numpy.floor(x))
            x2 = x1+1
            y1 = int(numpy.floor(y))
            y2 = y1+1
            V11 = cartesian_pattern[int(numpy.floor(y)),int(numpy.floor(x))]
            V12 = cartesian_pattern[int(numpy.floor(y)),int(numpy.floor(x))+1]
            V21 = cartesian_pattern[int(numpy.floor(y))+1,int(numpy.floor(x))]
            V22 = cartesian_pattern[int(numpy.floor(y))+1,int(numpy.floor(x))+1]
            polar_pattern[r,i_theta] = V11*(x2-x)*(y2-y) + V12*(x-x1)*(y2-y) + V21*(x2-x)*(y-y1) + V22*(x-x1)*(y-y1)
    return polar_pattern

def cartesian_to_radial(cartesian,N_theta):
    return numpy.mean(cartesian_to_polar(cartesian,N_theta),1)

def draw_circle(Nx,Ny,diameter):
    X,Y = numpy.meshgrid(numpy.arange(-Nx/2.0+0.5,Nx/2.0+0.5,1),numpy.arange(-Ny/2.0+0.5,Ny/2.0+0.5,1))
    circle = numpy.sqrt(X**2+Y**2)    
    circle[circle>diameter/2.0] = 0
    circle[circle!=0] = 1
    return circle

def downsample(array2d_raw,factor,mode="pick"):
    array2d = numpy.array(array2d_raw,dtype=array2d_raw.dtype)
    factor = int(factor)
    if factor == 1:
        return array2d
    available_modes = ["pick","integrate"]
    if not mode in available_modes:
        print "ERROR: %s is not a valid mode." % mode
        return 0
    Ny = array2d.shape[0]
    Nx = array2d.shape[1]
    if mode == "pick": 
        Ny_new = int(numpy.ceil(1.0*Ny/factor))
        Nx_new = int(numpy.ceil(1.0*Nx/factor))  
        array2d_new = numpy.zeros(Nx_new*Ny_new,dtype=array2d.dtype)  
        array2d_flat = array2d.flatten()
        for i in numpy.arange(0,Nx_new*Ny_new,1):
            ind = i%Nx_new*factor+(i/Nx_new)*Nx*factor
            array2d_new[i] = array2d_flat[ind]
        return numpy.reshape(array2d_new,(Ny_new,Nx_new))
    elif mode == "integrate":
        #Ny_new = int(numpy.floor(1.0*Ny/factor))
        #Nx_new = int(numpy.floor(1.0*Nx/factor))  
        Ny_new = int(numpy.ceil(1.0*Ny/factor))
        Nx_new = int(numpy.ceil(1.0*Nx/factor))  
        array2d_new = numpy.zeros(shape=(Ny_new,Nx_new),dtype=array2d.dtype)
        for y_new in numpy.arange(0,Ny_new,1):
            for x_new in numpy.arange(0,Nx_new,1):
                y_min = y_new*factor
                y_max = (y_new+1)*factor
                x_min = x_new*factor
                x_max = (x_new+1)*factor
                if y_max < Ny and x_max < Nx:
                    array2d_new[y_new,x_new] = array2d[y_min:y_max,x_min:x_max].mean()
                else:
                    if y_max >= Ny and x_max >= Nx:
                        array2d_new[y_new,x_new] = array2d[y_min:,x_min:].mean()
                    elif y_max >= Ny:
                        array2d_new[y_new,x_new] = array2d[y_min:,x_min:x_max].mean()
                    elif x_max >= Nx:
                        array2d_new[y_new,x_new] = array2d[y_min:y_max,x_min:].mean()
        return array2d_new

def crop(pattern,cropLength,center=None,bg=0,masking_threshold=None):
    if numpy.isscalar(cropLength):
        cropLength_x = cropLength
        cropLength_y = cropLength
    else:
        cropLength_x = cropLength[1]
        cropLength_y = cropLength[0]
    if not center:
        x_center = (pattern.shape[1] >> 1) - 0.5
        y_center = (pattern.shape[0] >> 1) - 0.5
    elif center == "center_of_mass":
        [y_center,x_center] = center_of_mass(pattern,masking_threshold)
        print [y_center,x_center]
    else:
        x_center = center[1]
        y_center = center[0]

    x_start = x_center-0.5-((cropLength_x >> 1)-1)
    x_stop = x_start + cropLength_x
    if x_start < 0:
        x_start = 0
        xc_start = -x_start
    else:
        xc_start = 0
    if x_stop > pattern.shape[1]:
        x_stop = pattern.shape[1] 
        xc_stop = x_stop - pattern.shape[1] 
    else:
        xc_stop = cropLength_x
    Nx = x_stop - x_start

    y_start = y_center-0.5-((cropLength_y >> 1)-1)
    y_stop = y_start + cropLength_y
    if y_start < 0:
        y_start = 0
        yc_start = -y_start
    else:
        yc_start = 0
    if y_stop > pattern.shape[0]:
        y_stop = pattern.shape[0] 
        yc_stop = y_stop - pattern.shape[0] 
    else:
        yc_stop = cropLength_y
    Ny = y_stop - y_start

    patternCropped = numpy.ones(shape=(cropLength_y,cropLength_x),dtype=pattern.dtype)*bg
    patternCropped[yc_start:yc_stop,xc_start:xc_stop] = pattern[y_start:y_stop,x_start:x_stop]
    return patternCropped

def turncw(array2d):
    array2d_turned = numpy.zeros_like(array2d)
    for x in range(0,len(array2d[0])):
        temp_list=list(array2d[:,x])
        temp_list.reverse()
        array2d_turned[x,:] = numpy.array(temp_list,dtype=array2d.dtype).T
    return array2d_turned

def turnccw(array2d):
    array2d_turned = numpy.zeros(shape=(array2d.shape[1],array2d.shape[0]),dtype=array2d.dtype)
    N = len(array2d_turned)-1
    for x in range(0,len(array2d[0])):
        array2d_turned[N-x,:] = array2d[:,x].T
    return array2d_turned

def horizontalmirr(array2d):
    array2d_mirrored = list(array2d.copy())
    array2d_mirrored.reverse()
    array2d_mirrored = numpy.array(array2d_mirrored)
    return array2d_mirrored

# only for Nx=Ny=Nz
#def slice(dm3d,phi,theta,psi):
#    #voxelcoordmatrix = get_voxel_coord_matrix(Nx,Ny,Nz)
#    N = dm3d.shape[0]
#    dm2dslice = numpy.zeros(N**2)
#    X2,X1 = numpy.meshgrid(numpy.arange(-(N-1)/2.0,(N/2-1)/2.0,1.0),numpy.arange(-(N-1)/2.0,(N/2-1)/2.0,1.0))
#    X0 = numpy.zeros_like(X1)
#    coord_slice = numpy.array([X0,X1,X2])
#    voxelcoordslice = _rotate_voxel_coord_slice(phi,theta,psi)
#    for i in range(0,len(dm2dslice)):
#        coord = voxelcoordslice[i]
#        dm2dslice[i] = interpolate3d(dm3d,coord)

#def interpolate3d(dm3d,coord,interpolation="linear"):
#    x0 = coord[0]
#    N0 = dm3d.shape[0]
#    x1 = coord[1]
#    N1 = dm3d.shape[1]
#    x2 = coord[2]
#    N2 = dm3d.shape[2]
#    if x0 > N0-1 or x1 > N1-1 or x2 > N2-1 or x0 < 0 or x1 < 0 or x2 < 0:
#        value = 0
#    else:
#        if interpolation == "linear":
#            value = 0
#            cases = [numpy.floor,lambda x: 1.0 + numpy.floor(x)]
#            for x0_func in cases:
#                for x1_func in cases:
#                    for x2_func in cases:
#                        value +=\
#                            (1.0-abs(x0_func(x0) - x0))*\
#                            (1.0-abs(x1_func(x1) - x1))*\
#                            (1.0-abs(x2_func(x2) - x2))*\
#                            dm3d[int(x0_func(x0)),int(x1_func(x1)),int(x2_func(x2))]
#        elif interpolation == "nn":
#            x0_rounded = numpy.round_(x0) 
#            x1_rounded = numpy.round_(x1) 
#            x2_rounded = numpy.round_(x2) 
#            value = dm3d[int(x0_rounded),int(x1_rounded),int(x2_rounded)]
#    return value

def pad_zeros(arr_orig,factor,shifted=False):
    arr = arr_orig.copy()
    Ny = arr.shape[0]
    Nx = arr.shape[1]
    Ny_new = int(round(Ny*factor)) 
    Nx_new = int(round(Nx*factor))
    Nx_d = Nx_new-Nx
    Ny_d = Ny_new-Ny
    arr_out = numpy.zeros(shape=(Ny_new,Nx_new),dtype=arr.dtype)
    if shifted:
        arr = numpy.fft.fftshift(arr)
    arr_out[Ny_d/2:Ny_d/2+Ny,Nx_d/2:Nx_d/2+Nx] = arr[:,:]
    if shifted:
        arr_out = numpy.fft.fftshift(arr_out)
    return arr_out
    
def depixelate(arr_orig,factor,shifted=False):
    arr = arr_orig.copy()
    farr = numpy.fft.fft2(arr)
    if shifted:
        pfarr = pad_zeros(farr,factor,False)
    else:
        pfarr = pad_zeros(farr,factor,True)
    parr = numpy.fft.ifft2(numpy.fft.fftshift(pfarr))
    return parr

def splitx(img,x_split=None):
    if x_split == None:
        x_split = img.shape[1]/2-0.5
    p1 = img[:,:x_split+0.5]
    p2 = img[:,x_split+0.5:]
    return [p1,p2]

def splity(img,y_split=None):
    if y_split == None:
        y_split = img.shape[0]/2-0.5
    p1 = img[:y_split+0.5,:]
    p2 = img[y_split+0.5:,:]
    return [p1,p2]

# Assume vertical slit (x = slitposition)
def stitch(img_in,splitposition,cx,cy,dx,dy):
    [p1,p2] = splitx(img_in,splitposition)
    Ny_out = img_in.shape[0]+abs(dy)
    Nx_out = img_in.shape[1]+abs(dx)
    Nx_p1 = p1.shape[1]
    Ny_p = p1.shape[0]
    img_out = numpy.zeros(shape=(Ny_out,Nx_out))
    if dy>=0:
        img_out[:Ny_p,:Nx_p1] = p1[:,:]
        img_out[dy:,Nx_p1+dx:] = p2[:,:]
    else:
        img_out[-dy:,:Nx_p1] = p1[:,:]
        img_out[:Ny_p,Nx_p1+dx:] = p2[:,:]
    return img_out

def resize2d(arr2d,dX_old,dX_new,X_new=None):
    from scipy.interpolate import interp2d
    from numpy import linspace
    if not X_new: X_new = arr2d.shape[1]
    x,y = numpy.meshgrid(numpy.arange(0,arr2d.shape[0])*dX_old,numpy.arange(0,arr2d.shape[1])*dX_old)
    newfunc = interp2d(x,y,arr2d,fill_value=0.0,kind='cubic')
    x_new = numpy.linspace(0,X_new*dX_new,dX_old/dX_new)
    y_new = numpy.linspace(0,X_new*dX_new,dX_old/dX_new)
    map2d_resized = newfunc(x_new,y_new)
    return map2d_resized

def interpolate3d(arr3d,factor):
    import enthought.mayavi.mlab as m
    N = arr3d.shape[0]
    arr3d_new = arr3d.copy()
    farr3d = numpy.fft.fftn(arr3d_new)
    farr3d = numpy.fft.fftshift(farr3d)
    farr3d_new = numpy.zeros(shape=(N*factor,N*factor,N*factor),dtype="complex")
    farr3d_new[round((N*factor-N)/2.0):round((N*factor-N)/2.0)+N,
               round((N*factor-N)/2.0):round((N*factor-N)/2.0)+N,
               round((N*factor-N)/2.0):round((N*factor-N)/2.0)+N] = farr3d[:,:,:]
    farr3d = farr3d_new
    farr3d = numpy.fft.fftshift(farr3d)
    arr3d_new = numpy.fft.ifftn(farr3d)
    return arr3d_new

def smooth3d(arr3d,factor):
    N = arr3d.shape[0]
    arr3d_new = arr3d.copy()
    farr3d = numpy.fft.fftn(arr3d_new)
    farr3d = numpy.fft.fftshift(farr3d)
    X,Y,Z = numpy.mgrid[-N/2:-N/2+N,-N/2:-N/2+N,-N/2:-N/2+N]
    Rsq = (X**2+Y**2+Z**2)
    kernel = numpy.exp(-Rsq/(N/2./(1.*factor))**2)
    farr3d = numpy.fft.fftshift(farr3d*kernel)
    arr3d_new = numpy.fft.ifftn(farr3d)
    return arr3d_new

def downsample3d_fourier(arr3d,factor):
    N = arr3d.shape[0]
    N_new = round(N*factor/2.0)*2
    arr3d_new = arr3d.copy()
    farr3d = numpy.fft.fftn(arr3d_new)
    farr3d = numpy.fft.fftshift(farr3d)
    A = farr3d.sum()
    farr3d = farr3d[(N-N*factor)/2:(N-N*factor)/2+N_new,(N-N*factor)/2:(N-N*factor)/2+N_new,(N-N*factor)/2:(N-N*factor)/2+N_new]
    B = farr3d.sum()
    farr3d /= (N/(1.0*N_new))**3.0
    farr3d = numpy.fft.fftshift(farr3d)
    arr3d_new = numpy.fft.ifftn(farr3d)
    return arr3d_new

def interpolate2d(arr2d,factor):
    N = arr2d.shape[0]
    arr2d_new = arr2d.copy()
    farr2d = numpy.fft.fftn(arr2d_new)
    farr2d = numpy.fft.fftshift(farr2d)
    farr2d_new = numpy.zeros(shape=(N*factor,N*factor),dtype="complex")
    farr2d_new[round((N*factor-N)/2.0):round((N*factor-N)/2.0)+N,
               round((N*factor-N)/2.0):round((N*factor-N)/2.0)+N] = farr2d[:,:]
    farr2d = farr2d_new
    farr2d = numpy.fft.fftshift(farr2d)
    arr2d_new = numpy.fft.ifftn(farr2d)
    return arr2d_new

def cut_edges(potential_cut_positions,normal_vectors,radius,dX,s=2.0):
    a = radius*(16*numpy.pi/5.0/(3+numpy.sqrt(5)))**(1/3.0)
    Rmax = numpy.sqrt(10.0+2*numpy.sqrt(5))*a/4.0 # radius at corners
    Rmin = numpy.sqrt(3)/12*(3.0+numpy.sqrt(5))*a # radius at faces
    nRmax = Rmax/dX
    nRmin = Rmin/dX
    N = int(numpy.ceil(2*(nRmax+1.0)))
    r_pix = dX*(3/(4*numpy.pi))**(1/3.0)
    #s = 2.0
    cutmap = numpy.ones(len(potential_cut_positions))
    for m in range(0,len(normal_vectors)):
        normal_vector = normal_vectors[m]
        for i in range(0,len(potential_cut_positions)):
            [iz,iy,ix] = potential_cut_positions[i]
            rsq = (iz-N/2.0-0.5)**2+(iy-N/2.0-0.5)**2+(ix-N/2.0-0.5)**2
            if rsq < (nRmin-s/2.0)**2:
                pass
            elif rsq > (nRmax+s/2.0)**2:
                cutmap[i] = 0.0
            elif cutmap[i] != 0.0:
                r = numpy.array([iz-N/2.0-0.5,iy-N/2.0-0.5,ix-N/2.0-0.5])
                delta = numpy.dot(1.0*r,1.0*normal_vector)/numpy.sqrt(numpy.dot(1.0*normal_vector,1.0*normal_vector)) - nRmin
                if delta > s/2.0:
                    cutmap[i] = 0.0
                elif abs(delta) <= s/2.0:
                    cutmap[i] = cutmap[i]*(0.5-delta/s)
    return cutmap

def get_random_circle_positions(N,d,dimension=2):
    # position circles randomly in a cube [-0.5..0.5,0.5..0.5,-0.5..0.5] without allowing overlap with border and neighboring circles (if d > 0)
    X = numpy.zeros(shape=(N,dimension))
    i = 0
    n_fails = 0
    dsq = d**2
    while i < N:
        X[i][:] = numpy.random.rand(dimension)[:]
        intersect = False
        for j in range(i):
            if numpy.dot(X[i]-X[j],X[i]-X[j]) <= dsq:
                intersect = True
                break
        if intersect == False:
            i += 1
        else:
            n_fails += 1
        if n_fails > 1E4:
            print "ERROR: Find no place to put more circles onto plane."
            return None
    return X


    
def rotate_3d_grid(X,Y,Z,eul_ang0,eul_ang1,eul_ang2):
    if eul_ang0 != 0.0 or eul_ang1 != 0.0 or eul_ang2 != 0.0:
        sizeX = X.shape[2]
        sizeY = X.shape[1]
        sizeZ = X.shape[0]   
        for xi in numpy.arange(0,sizeX,1.0):
            config.OUT.write("%i/%i\n" % (xi+1,sizeX))
            for yi in numpy.arange(0,sizeY,1.0):
                for zi in numpy.arange(0,sizeZ,1.0):
                    new_vector = pengtools.rotation(numpy.array([Z[zi,yi,xi],Y[zi,yi,xi],X[zi,yi,xi]]),eul_ang0,eul_ang1,eul_ang2)
                    X[zi,yi,xi] = new_vector[2]
                    Y[zi,yi,xi] = new_vector[1]
                    Z[zi,yi,xi] = new_vector[0]
    return [X,Y,Z]


def get_icosahedron_normal_vectors(euler_1=0.,euler_2=0.,euler_3=0.):
    # construct normal vectors of faces
    phi = (1+numpy.sqrt(5))/2.0
    ri = phi**2/2./numpy.sqrt(3.)
    # normal vectors for every vertice
    x1 = numpy.array([0.0,1.0,phi])
    x2 = numpy.array([0.0,1.0,-phi])
    x3 = numpy.array([0.0,-1.0,phi])
    x4 = numpy.array([0.0,-1.0,-phi]) 
    x5 = numpy.array([1.0,phi,0.0])
    x6 = numpy.array([1.0,-phi,0.0])
    x7 = numpy.array([-1.0,phi,0.0])
    x8 = numpy.array([-1.0,-phi,0.0])
    x9 = numpy.array([phi,0.0,1.0])
    x10 = numpy.array([-phi,0.0,1.0])
    x11 = numpy.array([phi,0.0,-1.0])
    x12 = numpy.array([-phi,0.0,-1.0])
    X = [x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12]
    # angle between normals
    an = round(numpy.dot(x5,x1))

    def cont_element(el,l):
        for i in range(0,len(l)):
            if (el == l[i]).all():
                return True
        return False

    def angles_match(y1,y2,y3):
        if round(numpy.dot(y1,y2)) == an and round(numpy.dot(y2,y3)) == an and round(numpy.dot(y3,y1)) == an:
            return True
        else:
            return False

    n_list = []
    for i in range(0,len(X)):
        for j in range(0,len(X)):
            for k in range(0,len(X)):
                n = (X[i]+X[j]+X[k])/6./ri
                if angles_match(X[i],X[j],X[k]) and not cont_element(n,n_list):
                    n_list.append(n)

                
    if euler_1 != 0. or euler_2 != 0. or euler_3 != 0.:
        for i in range(0,len(n_list)):
            n_list[i] = pengtools.rotation(n_list[i],euler_1,euler_2,euler_3)


    return n_list
