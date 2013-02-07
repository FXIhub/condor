import os,re,sys,h5py,pylab,numpy,time
import config
import proptools

def radial_pixel_average(image,**kargs):
    if 'cx' in kargs: cx = kargs['cx']
    else: cx = (image.shape[1]-1)/2.0
    if 'cy' in kargs: cy = kargs['cy'] 
    else: cy = (image.shape[0]-1)/2.0
    x = pylab.arange(0,image.shape[1],1.0)-cx
    y = pylab.arange(0,image.shape[1],1.0)-cy
    X,Y = pylab.meshgrid(x,y)
    R = pylab.sqrt(X**2+Y**2)
    R = R.round()
    R[image==pylab.Inf] = -1
    radii = pylab.arange(R.min(),R.max()+1,1)
    values = pylab.zeros_like(radii)
    for i in range(0,len(radii)):
        values[i] = image[R==radii[i]].mean()
    if 'rout' in kargs: return pylab.array([radii,values])
    else:return values
        
def radial_pixel_sum(image,**kargs):
    if 'cx' in kargs: cx = kargs['cx']
    else: cx = (image.shape[1]-1)/2.0
    if 'cy' in kargs: cy = kargs['cy'] 
    else: cy = (image.shape[0]-1)/2.0
    x = pylab.arange(0,image.shape[1],1.0)-cx
    y = pylab.arange(0,image.shape[1],1.0)-cy
    X,Y = pylab.meshgrid(x,y)
    R = pylab.sqrt(X**2+Y**2)
    R = R.round()
    radii = pylab.arange(R.min(),R.max()+1,1)
    values = pylab.zeros_like(radii)
    for i in range(0,len(radii)):
        values[i] = image[R==radii[i]].sum()
    if 'rout' in kargs: return pylab.array([radii,values])
    else:return values

def cartesian_to_polar(cartesian_pattern,N_theta,x_center=None,y_center=None):
    Nx = cartesian_pattern.shape[1]
    Ny = cartesian_pattern.shape[0]
    R = int(min([Nx,Ny])/2.0-1)
    polar_pattern = pylab.zeros(shape=(R,N_theta))
    if not x_center:
        x_center = Nx/2.0-0.5
    if not y_center:
        y_center = Ny/2.0-0.5
    for i_theta in range(0,N_theta):
        for r in range(0,R):
            theta = 2*pylab.pi*i_theta/float(N_theta)
            x = x_center + r * pylab.sin(theta)
            y = y_center + r * pylab.cos(theta)
            # bilinear interpolation
            x1 = int(pylab.floor(x))
            x2 = x1+1
            y1 = int(pylab.floor(y))
            y2 = y1+1
            V11 = cartesian_pattern[int(pylab.floor(y)),int(pylab.floor(x))]
            V12 = cartesian_pattern[int(pylab.floor(y)),int(pylab.floor(x))+1]
            V21 = cartesian_pattern[int(pylab.floor(y))+1,int(pylab.floor(x))]
            V22 = cartesian_pattern[int(pylab.floor(y))+1,int(pylab.floor(x))+1]
            polar_pattern[r,i_theta] = V11*(x2-x)*(y2-y) + V12*(x-x1)*(y2-y) + V21*(x2-x)*(y-y1) + V22*(x-x1)*(y-y1)
    return polar_pattern

def cartesian_to_radial(cartesian,N_theta):
    return pylab.mean(cartesian_to_polar(cartesian,N_theta),1)

def draw_circle(Nx,Ny,diameter):
    X,Y = pylab.meshgrid(pylab.arange(-Nx/2.0+0.5,Nx/2.0+0.5,1),pylab.arange(-Ny/2.0+0.5,Ny/2.0+0.5,1))
    circle = pylab.sqrt(X**2+Y**2)    
    circle[circle>diameter/2.0] = 0
    circle[circle!=0] = 1
    return circle

def downsample(array2d_raw,factor,mode="pick"):
    array2d = pylab.array(array2d_raw,dtype=array2d_raw.dtype)
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
        Ny_new = int(pylab.ceil(1.0*Ny/factor))
        Nx_new = int(pylab.ceil(1.0*Nx/factor))  
        array2d_new = pylab.zeros(Nx_new*Ny_new,dtype=array2d.dtype)  
        array2d_flat = array2d.flatten()
        for i in pylab.arange(0,Nx_new*Ny_new,1):
            ind = i%Nx_new*factor+(i/Nx_new)*Nx*factor
            array2d_new[i] = array2d_flat[ind]
        return pylab.reshape(array2d_new,(Ny_new,Nx_new))
    elif mode == "integrate":
        #Ny_new = int(pylab.floor(1.0*Ny/factor))
        #Nx_new = int(pylab.floor(1.0*Nx/factor))  
        Ny_new = int(pylab.ceil(1.0*Ny/factor))
        Nx_new = int(pylab.ceil(1.0*Nx/factor))  
        array2d_new = pylab.zeros(shape=(Ny_new,Nx_new),dtype=array2d.dtype)
        for y_new in pylab.arange(0,Ny_new,1):
            for x_new in pylab.arange(0,Nx_new,1):
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
    if pylab.isscalar(cropLength):
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

    patternCropped = pylab.ones(shape=(cropLength_y,cropLength_x),dtype=pattern.dtype)*bg
    patternCropped[yc_start:yc_stop,xc_start:xc_stop] = pattern[y_start:y_stop,x_start:x_stop]
    return patternCropped

def turncw(array2d):
    array2d_turned = pylab.zeros_like(array2d)
    for x in range(0,len(array2d[0])):
        temp_list=list(array2d[:,x])
        temp_list.reverse()
        array2d_turned[x,:] = pylab.array(temp_list,dtype=array2d.dtype).T
    return array2d_turned

def turnccw(array2d):
    array2d_turned = pylab.zeros(shape=(array2d.shape[1],array2d.shape[0]),dtype=array2d.dtype)
    N = len(array2d_turned)-1
    for x in range(0,len(array2d[0])):
        array2d_turned[N-x,:] = array2d[:,x].T
    return array2d_turned

def horizontalmirr(array2d):
    array2d_mirrored = list(array2d.copy())
    array2d_mirrored.reverse()
    array2d_mirrored = pylab.array(array2d_mirrored)
    return array2d_mirrored

# only for Nx=Ny=Nz
#def slice(dm3d,phi,theta,psi):
#    #voxelcoordmatrix = get_voxel_coord_matrix(Nx,Ny,Nz)
#    N = dm3d.shape[0]
#    dm2dslice = pylab.zeros(N**2)
#    X2,X1 = pylab.meshgrid(pylab.arange(-(N-1)/2.0,(N/2-1)/2.0,1.0),pylab.arange(-(N-1)/2.0,(N/2-1)/2.0,1.0))
#    X0 = pylab.zeros_like(X1)
#    coord_slice = pylab.array([X0,X1,X2])
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
#            cases = [pylab.floor,lambda x: 1.0 + pylab.floor(x)]
#            for x0_func in cases:
#                for x1_func in cases:
#                    for x2_func in cases:
#                        value +=\
#                            (1.0-abs(x0_func(x0) - x0))*\
#                            (1.0-abs(x1_func(x1) - x1))*\
#                            (1.0-abs(x2_func(x2) - x2))*\
#                            dm3d[int(x0_func(x0)),int(x1_func(x1)),int(x2_func(x2))]
#        elif interpolation == "nn":
#            x0_rounded = pylab.round_(x0) 
#            x1_rounded = pylab.round_(x1) 
#            x2_rounded = pylab.round_(x2) 
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
    arr_out = pylab.zeros(shape=(Ny_new,Nx_new),dtype=arr.dtype)
    if shifted:
        arr = pylab.fftshift(arr)
    arr_out[Ny_d/2:Ny_d/2+Ny,Nx_d/2:Nx_d/2+Nx] = arr[:,:]
    if shifted:
        arr_out = pylab.fftshift(arr_out)
    return arr_out
    
def depixelate(arr_orig,factor,shifted=False):
    arr = arr_orig.copy()
    farr = pylab.fft2(arr)
    if shifted:
        pfarr = pad_zeros(farr,factor,False)
    else:
        pfarr = pad_zeros(farr,factor,True)
    parr = pylab.ifft2(pylab.fftshift(pfarr))
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
    img_out = pylab.zeros(shape=(Ny_out,Nx_out))
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
    x,y = pylab.meshgrid(pylab.arange(0,arr2d.shape[0])*dX_old,pylab.arange(0,arr2d.shape[1])*dX_old)
    newfunc = interp2d(x,y,arr2d,fill_value=0.0,kind='cubic')
    x_new = pylab.linspace(0,X_new*dX_new,dX_old/dX_new)
    y_new = pylab.linspace(0,X_new*dX_new,dX_old/dX_new)
    map2d_resized = newfunc(x_new,y_new)
    return map2d_resized

def interpolate3d(arr3d,factor):
    import enthought.mayavi.mlab as m
    N = arr3d.shape[0]
    arr3d_new = arr3d.copy()
    farr3d = pylab.fftn(arr3d_new)
    farr3d = pylab.fftshift(farr3d)
    farr3d_new = pylab.zeros(shape=(N*factor,N*factor,N*factor),dtype="complex")
    farr3d_new[round((N*factor-N)/2.0):round((N*factor-N)/2.0)+N,
               round((N*factor-N)/2.0):round((N*factor-N)/2.0)+N,
               round((N*factor-N)/2.0):round((N*factor-N)/2.0)+N] = farr3d[:,:,:]
    farr3d = farr3d_new
    farr3d = pylab.fftshift(farr3d)
    arr3d_new = pylab.ifftn(farr3d)
    return arr3d_new

def smooth3d(arr3d,factor):
    N = arr3d.shape[0]
    arr3d_new = arr3d.copy()
    farr3d = pylab.fftn(arr3d_new)
    farr3d = pylab.fftshift(farr3d)
    X,Y,Z = pylab.mgrid[-N/2:-N/2+N,-N/2:-N/2+N,-N/2:-N/2+N]
    Rsq = (X**2+Y**2+Z**2)
    kernel = pylab.exp(-Rsq/(N/2./(1.*factor))**2)
    farr3d = pylab.fftshift(farr3d*kernel)
    arr3d_new = pylab.ifftn(farr3d)
    return arr3d_new

def downsample3d_fourier(arr3d,factor):
    N = arr3d.shape[0]
    N_new = round(N*factor/2.0)*2
    arr3d_new = arr3d.copy()
    farr3d = pylab.fftn(arr3d_new)
    farr3d = pylab.fftshift(farr3d)
    A = farr3d.sum()
    farr3d = farr3d[(N-N*factor)/2:(N-N*factor)/2+N_new,(N-N*factor)/2:(N-N*factor)/2+N_new,(N-N*factor)/2:(N-N*factor)/2+N_new]
    B = farr3d.sum()
    farr3d /= (N/(1.0*N_new))**3.0
    farr3d = pylab.fftshift(farr3d)
    arr3d_new = pylab.ifftn(farr3d)
    return arr3d_new

def interpolate2d(arr2d,factor):
    N = arr2d.shape[0]
    arr2d_new = arr2d.copy()
    farr2d = pylab.fftn(arr2d_new)
    farr2d = pylab.fftshift(farr2d)
    farr2d_new = pylab.zeros(shape=(N*factor,N*factor),dtype="complex")
    farr2d_new[round((N*factor-N)/2.0):round((N*factor-N)/2.0)+N,
               round((N*factor-N)/2.0):round((N*factor-N)/2.0)+N] = farr2d[:,:]
    farr2d = farr2d_new
    farr2d = pylab.fftshift(farr2d)
    arr2d_new = pylab.ifftn(farr2d)
    return arr2d_new

def cut_edges(potential_cut_positions,normal_vectors,radius,dX,s=2.0):
    a = radius*(16*pylab.pi/5.0/(3+pylab.sqrt(5)))**(1/3.0)
    Rmax = pylab.sqrt(10.0+2*pylab.sqrt(5))*a/4.0 # radius at corners
    Rmin = pylab.sqrt(3)/12*(3.0+pylab.sqrt(5))*a # radius at faces
    nRmax = Rmax/dX
    nRmin = Rmin/dX
    N = int(pylab.ceil(2*(nRmax+1.0)))
    r_pix = dX*(3/(4*pylab.pi))**(1/3.0)
    #s = 2.0
    cutmap = pylab.ones(len(potential_cut_positions))
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
                r = pylab.array([iz-N/2.0-0.5,iy-N/2.0-0.5,ix-N/2.0-0.5])
                delta = pylab.dot(1.0*r,1.0*normal_vector)/pylab.sqrt(pylab.dot(1.0*normal_vector,1.0*normal_vector)) - nRmin
                if delta > s/2.0:
                    cutmap[i] = 0.0
                elif abs(delta) <= s/2.0:
                    cutmap[i] = cutmap[i]*(0.5-delta/s)
    return cutmap

def get_random_circle_positions(N,d,dimension=2):
    # position circles randomly in a cube [-0.5..0.5,0.5..0.5,-0.5..0.5] without allowing overlap with border and neighboring circles (if d > 0)
    X = pylab.zeros(shape=(N,dimension))
    i = 0
    n_fails = 0
    dsq = d**2
    while i < N:
        X[i][:] = pylab.rand(dimension)[:]
        intersect = False
        for j in range(i):
            if pylab.dot(X[i]-X[j],X[i]-X[j]) <= dsq:
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

def generate_random_colloid_planes(edge_pix,diameter_pix,thickness_pix,Np,Nppp):
    # edge_pix: edge length of planes in pixel
    # diameter_pix: diameter of particles in pixel
    # thickness_pix: thickness of layer in pixel
    # Np: number of planes
    # Nppp: number of particles per plane
    box_hight = pylab.ceil(Np*thickness_pix+diameter_pix)
    box_width = pylab.ceil(edge_pix+diameter_pix)
    diameter = diameter_pix/(1.*box_width)
    config.OUT.write("Generate box of %i x %i x %i voxels.\n" % (box_width,box_width,box_hight))
    P = pylab.zeros((box_width,box_width,box_hight))
    X, Y, Z = numpy.mgrid[0:box_width,0:box_width,0:box_hight]
    positions = get_random_circle_positions(Nppp,diameter)
    y = positions[:,0]
    x = positions[:,1]
    x = x*(box_width-diameter_pix)+diameter_pix/2.
    y = y*(box_width-diameter_pix)+diameter_pix/2.
    for p in range(Np):
        config.OUT.write("Plane %i/%i" % (p+1,Np))
        for i in range(Nppp):
            R = pylab.sqrt((Y-y[i])**2+(X-x[i])**2+(Z-thickness_pix*(p+0.5))**2)
            P[R<diameter_pix/2.-1.0] = 1
            P[abs(R-diameter_pix/2.)<=1.0] = ((R[abs(R-diameter_pix/2.)<=1.0]-diameter_pix/2.0)+1.0)/2.0
    return P

#def position_spheres_in_icosahedron(d,N):
#    n_list = get_icosahedron_normal_vectors()
#    r = 0.5*

    # position circles randomly on a plane [0..1,0..1] without allowing overlap with border and neighboring circles (if d > 0)
#    positions = zeros((N,3))
#    i = 0
#    n_fails = 0
#    dsq = d**2
#    while i < N:
#        v = pylab.rand(3)
#        inicosahedron = True
#        for n in n_list:
#            if pylab.dot(v,n) > 1.0:
#                inicosahedron = False
#                break
#        if inicosahedron == True:
#            intersect = False
#            for j in range(i):
#                if (v[:]-positions[i,:])**2 <= dsq:
#                    intersect = True
#                    break
#            if intersect == False:
#                positions[i,:] = v.copy()[:]
#                i += 1
#        else:
#            n_fails += 1
#        if n_fails > 1E4:
#            print "ERROR: Find no place to put more circles onto plane."
#            return positions
#    
#    return positions
    
    
def rotate_3d_grid(X,Y,Z,eul_ang0,eul_ang1,eul_ang2):
    if eul_ang0 != 0.0 or eul_ang1 != 0.0 or eul_ang2 != 0.0:
        sizeX = X.shape[2]
        sizeY = X.shape[1]
        sizeZ = X.shape[0]   
        for xi in pylab.arange(0,sizeX,1.0):
            config.OUT.write("%i/%i\n" % (xi+1,sizeX))
            for yi in pylab.arange(0,sizeY,1.0):
                for zi in pylab.arange(0,sizeZ,1.0):
                    new_vector = proptools.rotation(pylab.array([Z[zi,yi,xi],Y[zi,yi,xi],X[zi,yi,xi]]),eul_ang0,eul_ang1,eul_ang2)
                    X[zi,yi,xi] = new_vector[2]
                    Y[zi,yi,xi] = new_vector[1]
                    Z[zi,yi,xi] = new_vector[0]
    return [X,Y,Z]

def get_icosahedrally_symmetrical_vectors(v,ico_n5=[1]):
    # construct normal vectors of faces
    phi = (1+sqrt(5))/2.0
    ri = phi**2/2./sqrt(3.)
    normalization_factor = 1/sqrt(phi**2+1.)
    # normal vectors for every vertice
    x1 = array([0.0,1.0,phi])*normalization_factor
    x2 = array([0.0,1.0,-phi])*normalization_factor
    x3 = array([0.0,-1.0,phi])*normalization_factor
    x4 = array([0.0,-1.0,-phi])*normalization_factor
    x5 = array([1.0,phi,0.0])*normalization_factor
    x6 = array([1.0,-phi,0.0])*normalization_factor
    x7 = array([-1.0,phi,0.0])*normalization_factor
    x8 = array([-1.0,-phi,0.0])*normalization_factor
    x9 = array([phi,0.0,1.0])*normalization_factor
    x10 = array([-phi,0.0,1.0])*normalization_factor
    x11 = array([phi,0.0,-1.0])*normalization_factor
    x12 = array([-phi,0.0,-1.0])*normalization_factor
    list5folds = [x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12]
    # angle between normals
    an = dot(x5,x1)

    def cont_element(el,l):
        for i in range(0,len(l)):
            if (sqrt(((el-l[i])**2)).sum()<0.01).any():
                return True
        return False

    def angles_match(y1,y2,y3):
        if abs(dot(y1,y2) - an) < 0.01 and abs(dot(y2,y3) - an) < 0.01 and abs(dot(y3,y1) - an) < 0.01:
            return True
        else:
            return False

    list3folds = []
    for i in range(0,len(list5folds)):
        for j in range(0,len(list5folds)):
            for k in range(0,len(list5folds)):
                n = (list5folds[i]+list5folds[j]+list5folds[k])/3.
                if angles_match(list5folds[i],list5folds[j],list5folds[k]) and not cont_element(n,list3folds):
                    list3folds.append(n)

    if euler_1 != 0. or euler_2 != 0. or euler_3 != 0.:
        for i in range(0,len(list3folds)):
            list3folds[i] = proptools.rotation(list3folds[i],euler_1,euler_2,euler_3)

    return list3folds

def get_icosahedron_normal_vectors(euler_1=0.,euler_2=0.,euler_3=0.):
    # construct normal vectors of faces
    phi = (1+pylab.sqrt(5))/2.0
    ri = phi**2/2./pylab.sqrt(3.)
    # normal vectors for every vertice
    x1 = pylab.array([0.0,1.0,phi])
    x2 = pylab.array([0.0,1.0,-phi])
    x3 = pylab.array([0.0,-1.0,phi])
    x4 = pylab.array([0.0,-1.0,-phi]) 
    x5 = pylab.array([1.0,phi,0.0])
    x6 = pylab.array([1.0,-phi,0.0])
    x7 = pylab.array([-1.0,phi,0.0])
    x8 = pylab.array([-1.0,-phi,0.0])
    x9 = pylab.array([phi,0.0,1.0])
    x10 = pylab.array([-phi,0.0,1.0])
    x11 = pylab.array([phi,0.0,-1.0])
    x12 = pylab.array([-phi,0.0,-1.0])
    X = [x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12]
    # angle between normals
    an = round(pylab.dot(x5,x1))

    def cont_element(el,l):
        for i in range(0,len(l)):
            if (el == l[i]).all():
                return True
        return False

    def angles_match(y1,y2,y3):
        if round(pylab.dot(y1,y2)) == an and round(pylab.dot(y2,y3)) == an and round(pylab.dot(y3,y1)) == an:
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
            n_list[i] = proptools.rotation(n_list[i],euler_1,euler_2,euler_3)


    return n_list
        
def lanczos_interp(data,N_new_,a=2.):
    
    N_new = int(round(N_new_))
    dim = len(data.shape)
    N = data.shape[0]

    fdata = pylab.fftn(data)
    sfdata = pylab.fftshift(fdata)

    if N_new%2 == 1:
        c = N/2-1
    else:
        c = N/2

    if dim == 1:
        sfdata_cropped = sfdata[c-N_new/2:c-N_new/2+N_new]
        X = (1.*pylab.arange(N_new)-N_new/2)/((N_new-1)/2.)
        lanczos_kernel = (pylab.sinc(X*a)*pylab.sinc(X))

    elif dim == 2:
        if data.shape[0] != data.shape[1]:
            print "ERROR: Only accept data with equal dimensions."
            return

        sfdata_cropped = sfdata[c-N_new/2:c-N_new/2+N_new,
                                c-N_new/2:c-N_new/2+N_new]
        X,Y = pylab.meshgrid((pylab.arange(N_new)-N_new/2)/((N_new-1)*2.),
                             (pylab.arange(N_new)-N_new/2)/((N_new-1)*2.))
        lanczos_kernel = pylab.sinc(X*a)*pylab.sinc(X)*pylab.sinc(Y*a)*pylab.sinc(Y)

    elif dim == 3:
        if data.shape[0] != data.shape[1] or data.shape[1] != data.shape[2]:
            print "ERROR: Only accept data with equal dimensions."
            return

        sfdata_cropped = sfdata[c-N_new/2:c-N_new/2+N_new,
                                c-N_new/2:c-N_new/2+N_new,
                                c-N_new/2:c-N_new/2+N_new]
        X,Y,Z = pylab.mgrid[:N_new,:N_new,:N_new]
        X = (X-N_new/2)/(2.*(N_new-1))
        Y = (Y-N_new/2)/(2.*(N_new-1))
        Z = (Z-N_new/2)/(2.*(N_new-1))
        lanczos_kernel = \
            pylab.sinc(X*a)*pylab.sinc(X)* \
            pylab.sinc(Y*a)*pylab.sinc(Y)* \
            pylab.sinc(Z*a)*pylab.sinc(Z)

    sfdata_cropped *= lanczos_kernel
    fdata_cropped = pylab.fftshift(sfdata_cropped)
    ffdata = pylab.ifftn(fdata_cropped)
    norm_factor = N_new/(1.*N)
    ffdata *= (N_new**dim/(1.*N)**dim)
    return ffdata

