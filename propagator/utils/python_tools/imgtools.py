#======================#
# Python tools - image #
#======================# 
#
# Author: Max Hantke
# Email: maxhantke@gmail.com

import os,re,sys,h5py,numpy,time,cmath
import logging
logger = logging.getLogger("imgtools")
import cxitools

this_folder = os.path.dirname(__file__)


def get_phase(x):
    return numpy.angle(x)

def get_R_and_Theta_map(Nx,Ny,cx=None,cy=None):
    if not cx:
        cx = (Nx-1)/2.0
    if not cy:
        cy = (Ny-1)/2.0
    x = numpy.arange(0,Nx,1.0)-cx
    y = numpy.arange(0,Ny,1.0)-cy
    X,Y = numpy.meshgrid(x,y)
    R = numpy.sqrt(X**2+Y**2)
    R = R.round()
    Theta = numpy.arctan(-Y/(X+numpy.finfo('float64').eps))
    Theta[X<0] += numpy.pi
    Theta += numpy.pi/2.0
    #numpy.imsave("Theta.png" , Theta)
    #numpy.imsave("X.png" , X)
    #numpy.imsave("Y.png" , Y)
    return [R,Theta]

def cone_pixel_average(image,N_theta,cx=None,cy=None):
    [R,Theta] = get_R_and_Theta_map(image.shape[1],image.shape[0],cx,cy)
    R[numpy.isfinite(image) == False] = -1
    radii = numpy.arange(R.min(),R.max()+1,1)
    if radii[0] == -1:
        radii = radii[1:]
    values = numpy.zeros(shape=(N_theta,len(radii)))
    for j in range(0,N_theta):
        theta_min = j/(1.0*N_theta)*2.0*numpy.pi
        theta_max = (j+1)/(1.0*N_theta)*2.0*numpy.pi
        theta_center = (theta_max+theta_min)/2.0
        theta_interval = theta_max-theta_min
        theta_image = image[abs(Theta-theta_center)<=theta_interval/2.0]
        theta_R = R[abs(Theta-theta_center)<=theta_interval/2.0]
        for i in range(0,len(radii)):
            temp = theta_image[theta_R==radii[i]].copy()
            temp = temp[numpy.isfinite(temp)]
            values[j,i] = temp.mean()
    return [radii,values]

#def radial_pixel_average(image,cx=None,cy=None):
#    [radii,values] = cone_pixel_average(image,1,cx,cy)
#    return [radii,values[0,:]]


def center_of_mass(img0,shifted=False):
    debug = False
    img = abs(img0)
    img = img/(1.*img.sum()+numpy.finfo("float32").eps)
    d = len(list(img.shape))
    cm = numpy.zeros(d)
    f = numpy.indices(img.shape)
    for i in range(d):
        f[i] = f[i]-numpy.ceil(img.shape[i]/2.)
        if not shifted:
            f[i,:] = numpy.fft.fftshift(f[i,:])[:]
        cm[i] = (f[i,:]*img[:]).sum()
        if debug:
            import pylab
            pylab.imsave(this_folder+"/testdata/f%i.png" % i,f[i])
    return cm

def cone_pixel_average_new(image,mask,N_theta,cx=None,cy=None,rdownsample=1):
    [R,Theta] = get_R_and_Theta_map(image.shape[1],image.shape[0],cx,cy)
    radii = numpy.arange(0,R.max()+1,1*rdownsample)
    values = numpy.zeros(shape=(N_theta,len(radii)))
    Mvalues = numpy.zeros(shape=(N_theta,len(radii)))
    for j in range(0,N_theta):
        theta_min = j/(1.0*N_theta)*2.0*numpy.pi
        theta_max = (j+1)/(1.0*N_theta)*2.0*numpy.pi
        theta_center = (theta_max+theta_min)/2.0
        theta_interval = theta_max-theta_min
        Mtheta =abs(Theta-theta_center)<=theta_interval/2.0
        Mtheta *= mask
        if Mtheta.sum() > 0:
            theta_image = image[Mtheta]
            theta_R = R[Mtheta]
            for i in range(0,len(radii)):
                m = abs(theta_R-radii[i])<=0.5*rdownsample
                if m.sum() > 0:
                    values[j,i] = theta_image[m].mean()
                    Mvalues[j,i] = 1
    return [values,Mvalues]

def radial_pixel_sum(image,cx=None,cy=None):
    if not cx:
        cx = (image.shape[1]-1)/2.0
    if not cy:
        cy = (image.shape[0]-1)/2.0
    x = numpy.arange(0,image.shape[1],1.0)-cx
    y = numpy.arange(0,image.shape[1],1.0)-cy
    X,Y = numpy.meshgrid(x,y)
    R = numpy.sqrt(X**2+Y**2)
    R = R.round()
    radii = numpy.arange(R.min(),R.max()+1,1)
    values = numpy.zeros_like(radii)
    for i in range(0,len(radii)):
        values[i] = image[R==radii[i]].sum()
    return numpy.array([radii,values])


def cartesian_to_polar(data,N_theta=100,cx0=None,cy0=None,r_max0=None,order=0):
    import scipy.ndimage

    if cx0 == None:
        cx = (data.shape[1]-1)/2.
    else:
        cx = cx0
    if cy0 == None:
        cy = (data.shape[0]-1)/2.
    else:
        cy = cy0

    if r_max0 == None:
        r_max = data.shape[0]/2
    else:
        r_max = r_max0

    def cartesian2polar(outcoords):
        ri,thetai = outcoords
        theta = thetai/(1.*(N_theta-1))*2*numpy.pi
        X = cx+ri*numpy.cos(theta)
        Y = cy+ri*numpy.sin(theta)
        return (Y,X)

    data_polar = scipy.ndimage.geometric_transform(data, cartesian2polar,order=order,output_shape=(r_max,N_theta))
    return data_polar


def cartesian_to_radial(cartesian,N_theta):
    return numpy.mean(cartesian_to_polar(cartesian,N_theta),1)

def draw_circle(Nx,Ny,diameter):
    X,Y = numpy.meshgrid(numpy.arange(-Nx/2.0+0.5,Nx/2.0+0.5,1),numpy.arange(-Ny/2.0+0.5,Ny/2.0+0.5,1))
    circle = numpy.sqrt(X**2+Y**2)    
    circle[circle>diameter/2.0] = 0
    circle[circle!=0] = 1
    return circle

def gaussian_smooth(I,sm,precision=1.):
    N = 2*int(numpy.round(precision*sm))+1
    if len(I.shape) == 2:
        import scipy.signal
        X,Y = numpy.meshgrid(numpy.arange(0,N,1),numpy.arange(0,N,1))
        X = X-N/2
        Y = Y-N/2
        R = numpy.sqrt(X**2 + Y**2)
        kernel = numpy.exp(R**2/(2.0*sm**2))
        kernel[abs(R)>N/2] = 0.
        kernel /= kernel.sum()
        Ism = scipy.signal.convolve2d(I,kernel,mode='same',boundary='fill')
        return Ism
    elif len(I.shape) == 1:
        X = numpy.arange(0,N,1)
        X = X-(N-1)/2.
        kernel = numpy.exp(X**2/(2.0*sm**2))
        kernel /= kernel.sum()
        Ism = numpy.convolve(I,kernel,mode='same')
        return Ism

def gaussian_smooth_2d1d(I,sm,precision=1.):
    N = 2*int(numpy.round(precision*sm))+1
    if len(I.shape) == 2:
        import scipy.signal
        kernel = numpy.zeros(shape=(N,N))
        X,Y = numpy.meshgrid(numpy.arange(0,N,1),numpy.arange(0,N,1))
        X = X-N/2
        kernel = numpy.exp(X**2/(2.0*sm**2))
        kernel /= kernel.sum()
        Ism = scipy.signal.convolve2d(I,kernel,mode='same',boundary='wrap')
        return Ism
    elif len(I.shape) == 1:
        print "Error input"
        return []

def gaussian_sharpening(image,sigma):
    imagefourier = numpy.fft2(image)
    Ny = image.shape[0]
    Nx = image.shape[1]
    X,Y = numpy.meshgrid(numpy.arange(-Nx/2.0+0.5,Nx/2.0+0.5,1.0),numpy.arange(-Ny/2.0+0.5,Ny/2.0+0.5,1.0))
    gauss = 1/numpy.sqrt(2*numpy.pi*sigma**2)*numpy.exp(-(X**2+Y**2)/(2*sigma**2))
    gauss = shift(gauss)
    imagefourier *= (1.0-gauss)
    return numpy.ifft2(imagefourier)
    
def downsample(array2d0,factor0,mode="pick",mask2d0=None,bad_bits=None,min_N_pixels=1):
    available_modes = ["pick","integrate"]#,"interpolate"]
    if not mode in available_modes:
        print "ERROR: %s is not a valid mode." % mode
        return
    factor = int(round(factor0))
    if factor == 1:
        if mask2d0 == None:
            return array2d0
        else:
            return [array2d0,mask2d0]
    array2d = numpy.array(array2d0,dtype=array2d0.dtype)
    if mask2d0 == None:
        mask2d = None
    else:
        mask2d = numpy.array(mask2d0,dtype="int16")
    Nx = array2d.shape[1]
    Ny = array2d.shape[0]
    if mode == "pick": 
        Y,X = numpy.indices(array2d0.shape)
        pick = ((Y%factor == 0)*(X%factor == 0))
        print pick.shape
        Ny_new = pick[:,0].sum()
        Nx_new = pick[0,:].sum()
        pick = pick.flatten()
        A = array2d.flatten()
        array2d_new = (A[pick]).reshape((Ny_new,Nx_new))
        if mask2d != None:
            M = mask2d.flatten()
            mask2d_new = (M[pick]).reshape((Ny_new,Nx_new))
            return [array2d_new,mask2d_new]
        else:
            return array2d_new
    elif mode == "integrate": # non-conservative if mask is given
        Nx_new = int(numpy.ceil(Nx/factor))
        Ny_new = int(numpy.ceil(Ny/factor))
        Nx = Nx_new * factor
        Ny = Ny_new * factor
        A = numpy.zeros(shape=(Ny,Nx),dtype=array2d.dtype)
        A[:array2d.shape[0],:array2d.shape[1]] = array2d[:,:]
        A = A.flat
        Y,X = numpy.indices((Ny,Nx))
        Y = Y.flatten()
        X = X.flatten()
        Y /= factor
        X /= factor
        superp = Y*Nx+X
        superp_order = superp.argsort()
        A = A[superp_order]
        A = A.reshape((Nx_new*Ny_new,factor*factor))
        if mask2d == None:
            B = A.sum(1)
            return B.reshape((Ny_new,Nx_new))
        if mask2d != None:
            AM = numpy.zeros(shape=(Ny,Nx),dtype="int16")
            AM[:mask2d.shape[0],:mask2d.shape[1]] = mask2d[:,:]
            AM = AM.flat
            AM = AM[superp_order]
            AM = AM.reshape((Nx_new*Ny_new,factor*factor))
            if bad_bits == None:
                B = (A*AM).sum(1)
                BN = AM.sum(1)
                BM = BN != 0
            else:
                B = (A*((AM & bad_bits) == 0)).sum(1)
                BN = ((AM & bad_bits) == 0).sum(1)
                BM = AM[:,0]
                for i in range(1,factor*factor):
                    BM |= AM[:,i]
                BM[BN >= min_N_pixels] = BM[BN >= min_N_pixels] & ~bad_bits
            B[BN >= min_N_pixels] = B[BN >= min_N_pixels] * factor/numpy.float64(BN[BN >= min_N_pixels])
            return [B.reshape((Ny_new,Nx_new)),BM.reshape((Ny_new,Nx_new))]


def downsample3d(array_raw,factor,mask):
    array_cp = array_raw.copy()
    factor = int(factor)
    if factor == 1:
        return array_cp
    Nz = array_cp.shape[0]
    Ny = array_cp.shape[1]
    Nx = array_cp.shape[2]
    Nz_new = int(numpy.ceil(1.0*Nz/factor))
    Ny_new = int(numpy.ceil(1.0*Ny/factor))
    Nx_new = int(numpy.ceil(1.0*Nx/factor))  
    array_new = numpy.zeros(shape=(Nz_new,Ny_new,Nx_new),dtype=array_cp.dtype)
    for z_new in numpy.arange(0,Nz_new,1):
        for y_new in numpy.arange(0,Ny_new,1):
            for x_new in numpy.arange(0,Nx_new,1):
                z_min = z_new*factor
                z_max = min([(z_new+1)*factor,Nz])
                y_min = y_new*factor
                y_max = min([(y_new+1)*factor,Ny])
                x_min = x_new*factor
                x_max = min([(x_new+1)*factor,Nx])
                array_new[z_new,y_new,x_new] = array_cp[z_min:z_max,y_min:y_max,x_min:x_max].sum()/(1.*mask[z_min:z_max,y_min:y_max,x_min:x_max].sum())
    return array_new

def downsample_spi(img,factor,mode="pick",ds_msk=True):
    import spimage
    img_array_new = downsample(img.image,factor,mode)
    img_new = spimage.sp_image_alloc(img_array_new.shape[1],img_array_new.shape[0],1)
    img_new.image[:,:] = img_array_new[:,:]
    if ds_msk:
        msk_array_new = downsample(img.mask,factor,mode)
        img_new.mask[:,:] = msk_array_new[:,:]
    else:
        img_new.mask[:,:] = 1
    return img_new

def crop(pattern,cropLength,center='middle',bg=0):
    if center == 'middle':
        x_center = (pattern.shape[1] - 1)/2.
        y_center = (pattern.shape[0] - 1)/2.
        temp = pattern.copy()
    else:
        if center == "center_of_mass":
            [y_center,x_center] = center_of_mass(pattern)
            x_center = numpy.ceil(pattern.shape[0]/2.) + x_center
            y_center = numpy.ceil(pattern.shape[0]/2.) + y_center
        else:
            x_center = center[1]
            y_center = center[0]
        temp = recenter(pattern,x_center,y_center)

    x_start = (pattern.shape[1]-cropLength)/2
    y_start = (pattern.shape[1]-cropLength)/2
    x_stop = x_start+cropLength
    y_stop = y_start+cropLength

    patternCropped = numpy.ones(shape=(cropLength,cropLength),dtype=pattern.dtype)*bg
    patternCropped = temp[y_start:y_stop,x_start:x_stop]
    return patternCropped

# crop around centre to maximal centrosymmetric size
def crop_max_around_center(array2d,c0,c1):
    d0 = array2d.shape[0] - 1 - c0
    d1 = array2d.shape[1] - 1 - c1
    N_new = min([c0,c1,d0,d1])*2+1
    start0 = c0-(N_new-1)/2.0
    stop0 = start0+N_new
    start1 = c1-(N_new-1)/2.0
    stop1 = start1+N_new
    array2d_new = array2d[start0:stop0,start1:stop1]
    return array2d_new

def diameter_extrema(image):
    N = 32
    polimage = cartesian_to_polar(image,N)
    radii = numpy.zeros(N)
    for i in range(0,N):
        radii[i] = polimage[:,i].argmin()
    diameters = radii[:N/2]+radii[N/2:]
    diameter_max = diameters.max()
    diameter_min = diameters.min()
    return [diameter_min,diameter_max]

    
def shift(pattern,center=None):
    Ny = pattern.shape[0]
    Nx = pattern.shape[1]
    if not center:
        center = [Ny/2,Nx/2]
    center = numpy.array([int(numpy.ceil(center[0])),int(numpy.ceil(center[1]))])
    patternShifted = numpy.zeros(shape=pattern.shape)
    patternShifted[Ny-center[0]:,Nx-center[1]:] = pattern[:center[0],:center[1]]
    patternShifted[Ny-center[0]:,:Nx-center[1]] = pattern[:center[0],center[1]:]
    patternShifted[:Ny-center[0],Nx-center[1]:] = pattern[center[0]:,:center[1]]
    patternShifted[:Ny-center[0],:Nx-center[1]] = pattern[center[0]:,center[1]:]
    return patternShifted

def shift_back(pattern,center):
    Ny = pattern.shape[0]
    Nx = pattern.shape[1]
    new_center = [0,0]
    new_center[0] = Ny-center[0]-1
    new_center[1] = Nx-center[1]-1
    return shift(pattern,new_center)

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

def fft_turn180(I):
    return numpy.fft.fftn(numpy.fft.fftn(I))/(1.*I.size)

def turn180(img,cx=None,cy=None):
    if cx == None:
        cx1 = (img.shape[0]-1)/2
    if cy == None:
        cy1 = (img.shape[0]-1)/2
    cx1 = round(cx*2)/2.
    cy1 = round(cy*2)/2.
    Nx1 = int(2*min([cx1,img.shape[1]-1-cx1]))+1
    Ny1 = int(2*min([cy1,img.shape[0]-1-cy1]))+1
    y_start = int(round(cy1-(Ny1-1)/2.))
    y_stop = int(round(cy1+(Ny1-1)/2.))+1
    x_start = int(round(cx1-(Nx1-1)/2.))
    x_stop = int(round(cx1+(Nx1-1)/2.))+1
    img_new = numpy.zeros(shape=(img.shape[0],img.shape[1]),dtype=img.dtype)
    #img_new = img.copy()
    img_new[y_start:y_stop,x_start:x_stop] = turnccw(turnccw(img[y_start:y_stop,x_start:x_stop]))
    return img_new
    
def test_turn180():
    outdir = "testout_trun180"
    os.system("mkdir %s" % outdir)
    os.system("rm %s/*" % outdir)
    A = get_test_image()
    import pylab
    pylab.imsave("%s/image.png" % outdir,A)
    B = turn180(A,A.shape[1]/3.,2*A.shape[0]/3.)
    pylab.imsave("%s/image_turned.png" % outdir,B)

def slide(img,dx=0,dy=0):
    if dx == 0 and dy == 0: imgout = img.copy()
    else:
        imgout=numpy.zeros_like(img)
        if dx > 0 and dy > 0: imgout[dy:,dx:] = img[:-dy,:-dx]
        elif dx < 0 and dy < 0: imgout[:dy,:dx] = img[-dy:,-dx:]
        elif dx < 0 and dy > 0: imgout[dy:,:dx] = img[:-dy,-dx:]
        elif dx > 0 and dy < 0: imgout[:dy,dx:] = img[-dy:,:-dx]
        elif dx > 0: imgout[:,dx:] = img[:,:-dx]
        elif dy > 0: imgout[dy:,:] = img[:-dy,:]
        elif dx < 0: imgout[:,:dx] = img[:,-dx:]
        elif dy < 0: imgout[:dy,:] = img[-dy:,:]
    return imgout
                                                
def turn180_(array2d,cx=None,cy=None):
    array2d_turned = turnccw(turnccw(array2d))
    dcx = cx-(array2d.shape[1]-1)/2.
    dcy = cy-(array2d.shape[0]-1)/2.
    array2d_turned = slide(array2d_turned,2*dcx,2*dcy)
    return array2d_turned

def horizontalmirr(array2d):
    array2d_mirrored = list(array2d.copy())
    array2d_mirrored.reverse()
    array2d_mirrored = numpy.array(array2d_mirrored)
    return array2d_mirrored

# only for Nx=Ny=Nz
def slice(dm3d,phi,theta,psi):
    #voxelcoordmatrix = get_voxel_coord_matrix(Nx,Ny,Nz)
    N = dm3d.shape[0]
    dm2dslice = numpy.zeros(N**2)
    X2,X1 = numpy.meshgrid(numpy.arange(-(N-1)/2.0,(N/2-1)/2.0,1.0),numpy.arange(-(N-1)/2.0,(N/2-1)/2.0,1.0))
    X0 = numpy.zeros_like(X1)
    coord_slice = numpy.array([X0,X1,X2])
    voxelcoordslice = _rotate_voxel_coord_slice(phi,theta,psi)
    for i in range(0,len(dm2dslice)):
        coord = voxelcoordslice[i]
        dm2dslice[i] = interpolate3d(dm3d,coord)
        
#def _get_voxel_coord_slice(coord_slice,phi,theta,psi):
     ### continue programming here


def interpolate3d(dm3d,coord,interpolation="linear"):
    x0 = coord[0]
    N0 = dm3d.shape[0]
    x1 = coord[1]
    N1 = dm3d.shape[1]
    x2 = coord[2]
    N2 = dm3d.shape[2]
    if x0 > N0-1 or x1 > N1-1 or x2 > N2-1 or x0 < 0 or x1 < 0 or x2 < 0:
        value = 0
    else:
        if interpolation == "linear":
            value = 0
            cases = [numpy.floor,lambda x: 1.0 + numpy.floor(x)]
            for x0_func in cases:
                for x1_func in cases:
                    for x2_func in cases:
                        value +=\
                            (1.0-abs(x0_func(x0) - x0))*\
                            (1.0-abs(x1_func(x1) - x1))*\
                            (1.0-abs(x2_func(x2) - x2))*\
                            dm3d[int(x0_func(x0)),int(x1_func(x1)),int(x2_func(x2))]
        elif interpolation == "nn":
            x0_rounded = numpy.round_(x0) 
            x1_rounded = numpy.round_(x1) 
            x2_rounded = numpy.round_(x2) 
            value = dm3d[int(x0_rounded),int(x1_rounded),int(x2_rounded)]
    return value

def pad_zeros(arr_orig,factor,shifted=False):
    arr = arr_orig.copy()
    Ny = arr.shape[0]
    Nx = arr.shape[1]
    Ny_new = int(round(Ny*factor)) 
    Nx_new = int(round(Nx*factor))
    Nx_d = Nx_new-Nx
    Ny_d = Ny_new-Ny
    arr_out = numpy.zeros(shape=(Ny_new,Nx_new),dtype=arr.dtype)
    if not shifted:
        arr = numpy.fftshift(arr)
    arr_out[Ny_d/2:Ny_d/2+Ny,Nx_d/2:Nx_d/2+Nx] = arr[:,:]
    if not shifted:
        arr_out = numpy.fftshift(arr_out)
    return arr_out

def fourier_upsample(arr_orig,factor,shifted=False):
    if not shifted:
        arr = arr_orig.copy()
    else:
        arr = numpy.fft.fftshift(arr_orig)
    farr = numpy.fft.fftn(arr)
    pfarr = pad_zeros(farr,factor)
    ups_arr = numpy.fft.ifftn(pfarr)
    return ups_arr

def depixelate(arr_orig,factor,shifted=False):
    arr = arr_orig.copy()
    farr = numpy.fft2(arr)
    if shifted:
        pfarr = pad_zeros(farr,factor,False)
    else:
        pfarr = pad_zeros(farr,factor,True)
    parr = numpy.ifft2(numpy.fftshift(pfarr))
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
    newfunc = interp2d(x,y,arr2d,fill_value=0.0,kind='linear')
    x_new = numpy.linspace(0,X_new*dX_new,dX_old/dX_new)
    y_new = numpy.linspace(0,X_new*dX_new,dX_old/dX_new)
    map2d_resized = newfunc(x_new,y_new)
    return map2d_resized

def interpolate3d(arr3d,factor):
    import enthought.mayavi.mlab as m
    N = arr3d.shape[0]
    arr3d_new = arr3d.copy()
    farr3d = numpy.fftn(arr3d_new)
    farr3d = numpy.fftshift(farr3d)
    farr3d_new = numpy.zeros(shape=(N*factor,N*factor,N*factor),dtype="complex")
    farr3d_new[round((N*factor-N)/2.0):round((N*factor-N)/2.0)+N,
               round((N*factor-N)/2.0):round((N*factor-N)/2.0)+N,
               round((N*factor-N)/2.0):round((N*factor-N)/2.0)+N] = farr3d[:,:,:]
    farr3d = farr3d_new
    farr3d = numpy.fftshift(farr3d)
    arr3d_new = numpy.ifftn(farr3d)
    return arr3d_new

def smooth1d(arr1d,sm):
    N = 1+2*sm
    x = numpy.arange(0,N,1) - sm
    kernel = numpy.zeros(N)
    kernel = numpy.exp(x**2/(1.0*sm**2))
    arr1d_new = numpy.convolve(arr1d,kernel,'same')
    return arr1d_new

def smooth3d(arr3d,factor):
    import enthought.mayavi.mlab as m
    N = arr3d.shape[0]
    arr3d_new = arr3d.copy()
    farr3d = numpy.fftn(arr3d_new)
    farr3d = numpy.fftshift(farr3d)
    X,Y,Z = numpy.mgrid[-N/2:-N/2+N,-N/2:-N/2+N,-N/2:-N/2+N]
    #R = numpy.sqrt(X**2+Y**2+Z**2)
    farr3d[abs(X)>N/2] = 0
    farr3d[abs(Y)>N/2] = 0
    farr3d[abs(Y)>N/2] = 0
    farr3d = numpy.fftshift(farr3d)
    arr3d_new = numpy.ifftn(farr3d)
    return arr3d_new

def downsample3d_fourier(arr3d,factor):
    import enthought.mayavi.mlab as m
    N = arr3d.shape[0]
    N_new = round(N*factor/2.0)*2
    arr3d_new = arr3d.copy()
    farr3d = numpy.fftn(arr3d_new)
    farr3d = numpy.fftshift(farr3d)
    A = farr3d.sum()
    farr3d = farr3d[(N-N*factor)/2:(N-N*factor)/2+N_new,(N-N*factor)/2:(N-N*factor)/2+N_new,(N-N*factor)/2:(N-N*factor)/2+N_new]
    B = farr3d.sum()
    farr3d /= (N/(1.0*N_new))**3.0
    farr3d = numpy.fftshift(farr3d)
    arr3d_new = numpy.ifftn(farr3d)
    return arr3d_new

def downsample2d_fourier(arr2d,factor):
    N = arr2d.shape[0]
    N_new = round(N*factor/2.0)*2
    arr2d_new = arr2d.copy()
    farr2d = numpy.fftn(arr2d_new)
    farr2d = numpy.fftshift(farr2d)
    A = farr2d.sum()
    farr2d = farr2d[(N-N*factor)/2:(N-N*factor)/2+N_new,(N-N*factor)/2:(N-N*factor)/2+N_new]
    B = farr2d.sum()
    farr2d /= (N/(1.0*N_new))**2.0
    farr2d = numpy.fftshift(farr2d)
    arr2d_new = numpy.ifftn(farr2d)
    return arr2d_new

def interpolate2d(arr2d,factor):
    N = arr2d.shape[0]
    arr2d_new = arr2d.copy()
    farr2d = numpy.fftn(arr2d_new)
    farr2d = numpy.fftshift(farr2d)
    farr2d_new = numpy.zeros(shape=(N*factor,N*factor),dtype="complex")
    farr2d_new[round((N*factor-N)/2.0):round((N*factor-N)/2.0)+N,
               round((N*factor-N)/2.0):round((N*factor-N)/2.0)+N] = farr2d[:,:]
    farr2d = farr2d_new
    farr2d = numpy.fftshift(farr2d)
    arr2d_new = numpy.ifftn(farr2d)
    #numpy.figure()
    #pylab.imshow(numpy.log10(farr2d.real))
    #pylab.imshow(arr2d_new.real)
    #pylab.show()
    return arr2d_new

def interpolate1d(arr1d,factor):
    N = len(arr1d)
    arr1d_new = arr1d.copy()
    farr1d = numpy.fftn(arr1d_new)
    farr1d = numpy.fftshift(farr1d)
    farr1d_new = numpy.zeros(N*factor,dtype="complex")
    farr1d_new[round((N*factor-N)/2.0):round((N*factor-N)/2.0)+N] = farr1d[:]
    farr1d = farr1d_new
    farr1d = numpy.fftshift(farr1d)
    arr1d_new = numpy.ifftn(farr1d)*factor
    #numpy.figure()
    #pylab.imshow(pylab.log10(farr2d.real))
    #pylab.imshow(arr2d_new.real)
    #pylab.show()
    return arr1d_new

def put_besides(img1,img2):
    img3 = numpy.zeros(shape=(max([img1.shape[0],img2.shape[0]]),img1.shape[1]+img2.shape[1]))
    img3[:img1.shape[0],:img1.shape[1]] = img1[:,:]
    img3[:img2.shape[0]:,img1.shape[1]:] = img2[:,:]
    return img3

def _radial(image,mode="mean",**kwargs):
    if mode == "mean": f = numpy.mean
    elif mode == "sum": f = numpy.sum
    elif mode == "std": f = numpy.std
    elif mode == "median": f = numpy.median
    else:
        print "ERROR: No valid mode given for radial projection."
        return
    if 'cx' in kwargs: cx = kwargs['cx']
    else: cx = (image.shape[1]-1)/2.0
    if 'cy' in kwargs: cy = kwargs['cy'] 
    else: cy = (image.shape[0]-1)/2.0
    R = get_R_and_Theta_map(image.shape[1],image.shape[0],cx,cy)[0]
    R = R.round()
    R[numpy.isfinite(image)==False] = -1
    radii = numpy.arange(R.min(),R.max()+1,1)
    if radii[0] == -1:
        radii = radii[1:]
    values = numpy.zeros_like(radii)
    for i in range(0,len(radii)):
        values[i] = f(image[R==radii[i]])
    if 'rout' in kwargs: return numpy.array([radii,values])
    else: return values
def radial_sum(image,**kwargs):
    return _radial(image,"sum",**kwargs)
def radial_std(image,**kwargs):
    return _radial(image,"std",**kwargs)
def radial_mean(image,**kwargs):
    return _radial(image,"mean",**kwargs)
def radial_median(image,**kwargs):
    return _radial(image,"median",**kwargs)

def _radial_fast(image,mode="mean",**kwargs):
    if mode != "mean" and mode != "sum":
        print "ERROR: No valid mode given for radial projection."
        return
    if 'cx' in kwargs: cx = int(round(kwargs['cx']))
    else: cx = (image.shape[1]-1)/2
    if 'cy' in kwargs: cy = int(round(kwargs['cy']))
    else: cy = (image.shape[0]-1)/2
    X, Y = numpy.indices(image.shape) 
    R = numpy.sqrt((X-cx)**2+(Y-cy)**2)
    ind = numpy.argsort(R.flat)
    # radii-sorted data
    sRi = R.flat[ind].astype(numpy.int16)   
    simage = numpy.float64(image.flat[ind])
    # array that is 1 for all indices where radius changes in relation to previous index
    dR = (sRi[1:] - sRi[:-1]) != 0
    csimage = numpy.cumsum(simage)
    rsimage = numpy.ones(dR.sum()+1,dtype="float")
    rsimage[0] = csimage[0]
    rsimage[1:] = (csimage[1:])[dR] - (csimage[:-1])[dR]
    if mode == "sum":
        return rsimage
    elif mode == "mean":
        csones = numpy.cumsum(numpy.ones(len(csimage),dtype="float"))
        rN = numpy.zeros(dR.sum()+1,dtype="float")
        rN[0] = 0.
        rN[1:] = (csones[1:])[dR] - (csones[:-1])[dR]
        return rsimage/csones

def radial_sum_fast(image,**kwargs):
    return _radial_fast(image,"sum",**kwargs)
def radial_mean_fast(image,**kwargs):
    return _radial_fast(image,"mean",**kwargs)

#def recenter(I,cx,cy):
#    dx = int(pylab.ceil(cx-(I.shape[1]-1)/2.))
#    dy = int(pylab.ceil(cy-(I.shape[0]-1)/2.))
#    I_recentered = pylab.zeros_like(I)
#    if   dx > 0 and dy > 0:
#        I_recentered[:-dy,:-dx] = I[dy:,dx:]
#    elif dx > 0 and dy < 0:
#        I_recentered[-dy:,:-dx] = I[:dy,dx:]
#    elif dx < 0 and dy > 0:
#        I_recentered[:-dy,-dx:] = I[dy:,:dx]
#    elif dx < 0 and dy < 0:
#        I_recentered[-dy:,-dx:] = I[:dy,:dx]
#    elif dx == 0 and dy < 0:
#        I_recentered[-dy:,:] = I[:dy,:]
#    elif dx == 0 and dy > 0:
#        I_recentered[:-dy,:] = I[dy:,:]
#    elif dx > 0 and dy == 0:
#        I_recentered[:,:-dx] = I[:,dx:]
#    elif dx < 0 and dy == 0:
#        I_recentered[:,-dx:] = I[:,:dx]
#    elif dx == 0 and dy == 0:
#        I_recentered[:,:] = I[:,:]
#    return I_recentered

def recenter(I,cx,cy,order=1):
    dx = int(numpy.ceil(cx-(I.shape[1]-1)/2.))
    dy = int(numpy.ceil(cy-(I.shape[0]-1)/2.))
    return pixel_translation(I,[dy,dx],order)

def downsample_position(position,downsampling):
    return (position-(downsampling-1)/2.)/(1.*downsampling)

def pair_correlation(I,M0=None):
    if M0 == None:
        M = numpy.ones(shape=I.shape,dtype="bool")
    else:
        M = M0

    Nx = I.shape[1]
    Ny = I.shape[0]

    CI = numpy.zeros_like(I)
    CN = numpy.zeros_like(M)
    Im = numpy.zeros_like(I)
    
    for dy in range(Ny):
        if M[dy,:].sum() > 1:
            Iy = I[dy,:][M[dy,:]]
            Im[dy,:] = (I[dy,:]-Iy.mean())/(Iy.std()+numpy.finfo('float64').eps)
    for dx in range(Nx):
        for ix in range(Nx):
            m = M[:,ix]*M[:,(ix+dx)%Nx]
            CN[:,dx] += m[:]
            CI[:,dx] += Im[:,ix]*Im[:,(ix+dx)%Nx]*m[:]
    
    CM = CN != 0
    return [CI,CM]

def test_pair_correlation():
    if True:
        Nx = 101
        Ny = 111
        X,Y = numpy.meshgrid(numpy.arange(Nx),numpy.arange(Ny))
        img = numpy.sin(2*numpy.pi/(10.)*X)
        msk = numpy.ones(shape=img.shape,dtype="bool")
        CI,CM = pair_correlation(img,msk)
    import pylab
    pylab.clf()
    outfolder = "./testdata/"
    pylab.plot(CI.sum(0))
    pylab.savefig('%s/test_Cpolar.png' % outfolder)
    pylab.imsave('%s/C.png' % outfolder,CI*numpy.log10(10*CM))
    pylab.imsave('%s/img.png' % outfolder,img)
    pylab.imsave('%s/msk.png' % outfolder,msk)

def radial_pair_correlation(I,M0=None,N_theta=None,cx=None,cy=None):
    Ip = cartesian_to_polar(I,N_theta,cx,cy)
    if M0 != None:
        Mp = cartesian_to_polar(M0,N_theta,cx,cy) == 1
    else:
        Mp = None
    return pair_correlation(Ip,Mp)

def test_radial_pair_correlation():
    name = "pair_correlation"
    if False:
        Nx = 120
        Ny = 100
        cx = 40
        cy = 50
        R,Theta = get_R_and_Theta_map(Nx,Ny,cx,cy)
        img = R*numpy.sin(Theta*10.)
        msk = numpy.ones(shape=img.shape,dtype="bool")
        msk[30:50,:] = False
        CI,CM = radial_pair_correlation(img,msk,100,cx,cy)
    if True:
        import propagator as p
        I = p.Input()
        ds = 2
        I.detector.init_mask(Nx=1024,
                             Ny=1024,
                             y_gap_size_in_pixel=23,
                             x_gap_size_in_pixel=0,
                             hole_diameter_in_pixel=70)
        I.propagation.rs_oversampling = 2.
        I.set_sample_icosahedral_virus_map(225E-09)
        I.sample.set_random_orientation()
        O = p.propagator(I)
        img = numpy.poisson(O.get_intensity_pattern())
        msk = I.detector.mask
        cx =  I.detector.cx
        cy = I.detector.cy
        CI,CM = radial_pair_correlation(img,msk,cx,cy)
    Ip = cartesian_to_polar(img,100,cx,cy)
    import pylab
    pylab.imsave("testdata/%s_img.png" % name,img*numpy.log10(msk*10))
    pylab.imsave("testdata/%s_Ip.png" % name,Ip)
    pylab.imsave("testdata/%s_PC.png" % name,CI*numpy.log10(CM*10))
    pylab.imsave("testdata/%s_limg.png" % name,numpy.log10(img)*numpy.log10(msk*10))
    pylab.imsave("testdata/%s_lIp.png" % name,numpy.log10(Ip))
    pylab.imsave("testdata/%s_lPC.png" % name,numpy.log10(CI)*numpy.log10(CM*10))

    

downsample_position = lambda x,N,binsize: (x-(binsize-1)/2.)*(N/(1.*binsize)-1)/(1.*(N-binsize))
upsample_position = lambda x,N,binsize: x*(N*binsize-binsize)/(1.*(N-1))+(binsize-1)/2.

def get_test_image():
    import Image,os
    filename = os.path.dirname(os.path.realpath(__file__)) + "/testdata/testmax_gray.png"
    I = Image.open(filename)
    Nx,Ny = I.size
    D = numpy.array(I.getdata())[:]
    D=D.reshape((Ny,Nx))
    return D    

# NOT TESTED
def phase_diff(imgA,imgB):
    A = imgA.copy()
    B = imgB.copy()
    A = A%(2*numpy.pi)
    B = B%(2*numpy.pi)
    return A-B


# should be done with stsci.image package in the future
#def pixel_translation(A,t,method="linear",fill_value=0):
#    from scipy.interpolate import griddata
#    d = len(list(A.shape))
#    g = numpy.indices(A.shape)
#    gt = numpy.indices(A.shape)
#    gt = list(gt)
#    g = list(g)
#    for i in range(d):
#        gt[i] = ((gt[i]+t[i]) % A.shape[i])
#        g[i] = g[i].flatten()
#    gt = tuple(gt)
#    g = tuple(g)
#    return griddata(g,A.flatten(),gt,method=method,fill_value=fill_value)

def pixel_translation(A,t,order=1):
    if A.dtype == "complex64":
        return (1.*pixel_translation(A.real,t,order)+1.j*pixel_translation(A.imag,t,order))
    from scipy import ndimage
    d = len(list(A.shape))
    g = numpy.indices(A.shape)
    gt = numpy.indices(A.shape)
    for i in range(d):
        gt[i] = (gt[i]+t[i]) % A.shape[i]
    return ndimage.map_coordinates(A, gt, order=order)


# should be done with stsci.image package in the future
def upsample(A,f0,order=1):
    if A.dtype == "complex64":
        return (1.*upsample(A.real,f,order)+1.j*upsample(A.imag,f,order))
    from scipy import ndimage
    d = len(list(A.shape))
    new_shape = tuple((numpy.array(A.shape)*f0).round())
    gt = numpy.float64(numpy.indices(new_shape))
    for i in range(d):
        f = (new_shape[i]-1)/float(A.shape[i]-1)
        gt[i] = gt[i]/f
    return ndimage.map_coordinates(A, gt, order=order)


# t = [t0,t1,...] (transferred into python from libspimage)
def fourier_translation(A,t,rotation=False):
    fA = numpy.fft.fftn(A)
    d = len(list(fA.shape))
    f = numpy.indices(fA.shape)
    for i in range(d):
        f[i] = f[i]-numpy.ceil(fA.shape[i]/2.)
        f[i,:] = numpy.fft.fftshift(f[i,:])[:]
    tmp = 0
    for i,ti,fi in zip(range(d),t,f):
        tmp = tmp + 2*numpy.pi*fi[:,:]*ti/fA.shape[i]
    A_translated = numpy.fft.ifftn(fA*numpy.exp(-1.j*tmp))
    #print "%e" % (abs(A).sum()/abs(A_translated).sum())
    return A_translated

def fourier_translation_test():
    import pylab
    A = get_test_image()
    pylab.imsave("testdata/fourier_translation_test_A.png",A,cmap=pylab.cm.gray)
    B = fourier_translation(A,[45,34])
    pylab.imsave("testdata/fourier_translation_test_B.png",B,cmap=pylab.cm.gray,vmin=A.min(),vmax=B.max())

def recover_translation(imgA,imgB,enantio=False):
    debug = False
    imgfA = numpy.fft.fftn(imgA)
    imgfB = numpy.fft.fftn(imgB)
    d = len(list(imgfA.shape))
    f = numpy.indices(imgfA.shape)
    for i in range(d):
        f[i] = f[i]-numpy.ceil(imgfA.shape[i]/2.)
        f[i,:] = numpy.fft.fftshift(f[i,:])[:]
    if enantio == False:
        # Check superposition with image
        imgB_new = imgB
        c = abs(numpy.fft.ifftn(imgfA*imgfB.conj()))
        turned = False
    else:
        imgB_turned = fft_turn180(imgB)
        imgfB_turned = numpy.fft.fftn(imgB_turned)
        # Check superposition with normal and rotated image
        cc = [abs(numpy.fft.ifftn(imgfA*imgfB.conj())),abs(numpy.fft.ifftn(imgfA*imgfB_turned.conj()))]
        if debug:
            os.system("rm %s/testdata/RT*" % this_folder)
            pylab.imsave(this_folder+"/testdata/RT_imgB_turned.png",abs(numpy.fft.fftshift(imgB_turned)),vmin=0.,vmax=abs(imgB).max())
            pylab.imsave(this_folder+"/testdata/RT_imgB.png",abs(numpy.fft.fftshift(imgB)),vmin=0.,vmax=abs(imgB).max())
            pylab.imsave(this_folder+"/testdata/RT_imgA.png",abs(numpy.fft.fftshift(imgA)),vmin=0.,vmax=abs(imgB).max())
            pylab.imsave(this_folder+"/testdata/RT_imgfB_turned.png",abs(numpy.fft.fftshift(imgfB_turned)),vmin=0.,vmax=abs(imgfB).max())
            pylab.imsave(this_folder+"/testdata/RT_imgfB.png",abs(numpy.fft.fftshift(imgfB)),vmin=0.,vmax=abs(imgfB).max())
            pylab.imsave(this_folder+"/testdata/RT_imgfA.png",abs(numpy.fft.fftshift(imgfA)),vmin=0.,vmax=abs(imgfB).max())
            pylab.imsave(this_folder+"/testdata/RT_CC0_%i.png" % k,cc[0])
            pylab.imsave(this_folder+"/testdata/RT_CC1_%i.png" % k,cc[1])
        Mcc = numpy.array([cc[0].max(),cc[1].max()])
        i_max = Mcc.argmax()
        if i_max == 0:
            turned = False
            c = cc[0]
        else:
            turned = True
            c = cc[1]
    index_max = c.argmax()
    translation = []
    for i in range(d):
        translation.append(f[i,:].flatten()[index_max])
    translation = numpy.array(translation)
    return [translation,turned]

# This functions translates image b so that it's phases 
# are as close as possible to a.
# The translation is done in fourier space and both images
# should be in real space
# (transferred into python from libspimage)
def maximize_overlap(imgA0,imgB0,enantio=False):
    imgA = imgA0.copy()
    imgB = imgB0.copy()
    [translation,turned] = recover_translation(imgA,imgB,enantio)
    if turned: imgB = fft_turn180(imgB)
    imgB = fourier_translation(imgB,translation)
    return [imgB,translation,turned]

def maximize_overlap_test():
    A = get_test_image()
    A[:A.shape[0]/3,:] = 0
    A[2*A.shape[0]/3:,:] = 0
    A[:,:A.shape[1]/3] = 0
    A[:,2*A.shape[1]/3:] = 0
    A = A[:,:]
    t = [-43,-23]
    B0 = abs(fourier_translation(A,t))
    for i,B in zip(range(2),[B0,fft_turn180(B0)]):
        C = maximize_overlap(A,B,True)
        print "Difference %i: %f" % (i,abs(A-C).sum())
        pylab.imsave("testdata/maximize_overlap_test_%i_A.png" % i,A,cmap=pylab.cm.gray)
        pylab.imsave("testdata/maximize_overlap_test_%i_B.png" % i,B,cmap=pylab.cm.gray,vmin=A.min(),vmax=A.max())
        pylab.imsave("testdata/maximize_overlap_test_%i_AB.png" % i,C,cmap=pylab.cm.gray,vmin=A.min(),vmax=A.max())

# under construction
# maximize phase match in fourier space
def maximize_phase_match(imgA,imgB,enantio=False):
    from scipy.optimize import leastsq
    imgfA = numpy.fft.fftn(imgA)
    imgfB = numpy.fft.fftn(imgB)
    pimgfA = numpy.angle(imgfA)
    pimgfB = numpy.angle(imgfB)

    pdiff = pimgfA-pimgfB

    d = len(list(imgfA.shape))
    f = numpy.indices(imgfA.shape)
    for i in range(d):
        f[i] = f[i]-numpy.ceil(imgfA.shape[i]/2.)
        f[i,:] = numpy.fft.fftshift(f[i,:])[:]

    if d == 1:
        p = lambda v: v[0] + v[1]*f
    elif d ==2:
        p = lambda v: v[0] + v[1]*f[0,:,:] + v[2]*f[1,:,:]
    elif d ==3:
        p = lambda v: v[0] + v[1]*f[0,:,:] + v[2]*f[1,:,:] + v[3]*f[2,:,:]
    v0 = numpy.array([0.]*(d+1))
    err = lambda v: ((p(v)-pdiff)**2).sum()
    v1, success = leastsq(lambda v: ones(len(v))*err(v),v0)
    
    if enantio:
        imgfB_e = imgfB.conj()
        pimgfB_e = numpy.angle(pimgfB)
        err_e = lambda v: (abs(p(v)-pdiff_e)).sum()
        v1_e, success_e = leastsq(lambda v: ones(len(v))*err_e(v),v0)
        
        if err_e(v1_e) < err(v1):
            v1 = v1_res
            pimgfB_res = -pimgfB

    fres = abs(imgB)*numpy.exp(1.j*(pimgfB+p(v1)))
    res = numpy.fft.ifftn(fres)
    return [ref,fres]

# in fourier space subtract phase ramp obtained by leastsq, corresponding to translation in real space
# input: real space
# output: real space and fourier space data
def minimize_phase_ramp(img,shifted=False,periodic_boundary=False):
    from scipy.optimize import leastsq,fmin_cobyla
    debug = True
    imgf = numpy.fft.fftn(numpy.fft.fftshift(img))
    pimgf = numpy.angle(imgf)+numpy.pi

    d = len(list(imgf.shape))
    f = numpy.indices(imgf.shape)
    for i in range(d):
        f[i] = f[i]-numpy.ceil(imgf.shape[i]/2.)
        f[i,:] = numpy.fft.fftshift(f[i,:])[:]
    if d == 1:
        p = lambda v: (v[0] + v[1]*f) % (2*numpy.pi)
        v00 = pimgf[f[0]==0][0]
        v01 = pimgf[f[0]==1][0]-v00
        v0 = [v00,v01]
    elif d ==2:
        p = lambda v: (v[0] + v[1]*f[0,:,:] + v[2]*f[1,:,:])  % (2*numpy.pi)
        v00 = pimgf[(f[0]==0)*(f[1]==0)][0]
        v01 = pimgf[(f[0]==1)*(f[1]==0)][0]-v00
        v02 = pimgf[(f[0]==0)*(f[1]==1)][0]-v00
        v0 = [v00,v01,v02]
    elif d ==3:
        p = lambda v: v[0] + (v[1]*f[0,:,:] + v[2]*f[1,:,:] + v[3]*f[2,:,:])  % (2*numpy.pi)
        v00 = pimgf[(f[0]==0)*(f[1]==0)*(f[2]==0)][0]
        v01 = pimgf[(f[0]==1)*(f[1]==0)*(f[2]==0)][0]-v00
        v02 = pimgf[(f[0]==0)*(f[1]==1)*(f[2]==0)][0]-v00
        v03 = pimgf[(f[0]==0)*(f[1]==0)*(f[2]==1)][0]-v00
        v0 = [v00,v01,v02,v03]
    err = lambda v: ((pimgf-p(v))**2).sum()
    v1, success = leastsq(lambda v: numpy.ones(len(v))*err(v),v0)

    v2 = v1.copy()
    if periodic_boundary:
        for i in range(d):
            m1 = v1[i+1]*imgf.shape[i]/(2.*numpy.pi)
            m2 = round(v1[i+1]*imgf.shape[i]/(2.*numpy.pi))
            v2[i+1] = m2/m1*v1[i+1]

    resf = abs(imgf)*numpy.exp(1.j*(pimgf-p(v2)))
    res = numpy.fft.ifftn(resf)

    if shifted:
        resf = numpy.fft.fftshift(resf)
        res = numpy.fft.fftshift(res)

    translation = (v2[1:]-numpy.pi)/(2*numpy.pi)*numpy.array(pimgf.shape)
    if debug:
        pylab.imsave(this_folder+"/testdata/subtract_phase_ramp_img0.png",abs(img))
        pylab.imsave(this_folder+"/testdata/subtract_phase_ramp_img1.png",abs(res))
        pylab.imsave(this_folder+"/testdata/subtract_phase_ramp_img0t.png",abs(fourier_translation(img,translation)))
    return [res,translation]



# Minimize the difference between the phases of a and b by adding a constant phase to b.
# (transferred into python from libspimage)
def phase_match(imgA,imgB,weights=None): # typically weights = (abs(imgA)*abs(imgB))
    diff = numpy.angle(imgA)-numpy.angle(imgB)
    if weights == None:
        w = 1/(1.*len(diff.flatten()) + numpy.finfo('float64').eps)
    else:
        w = weights / (weights.sum() + numpy.finfo('float64').eps)
    return (diff*w).sum()


def prtf(imgs0,msks0,**kwargs):
    debug = False
    logger0 = kwargs.get("logger",None)
    K = numpy.random.randint(1000)
    enantio = kwargs.get("enantio",False)
    shifted = kwargs.get("shifted",True)
    center_result = kwargs.get("center_result",None)
    pixels_to_exclude = kwargs.get("pixels_to_exclude",None)
    do_maximize_overlap = kwargs.get("maximize_overlap",True)
    do_minimize_phase_ramp = kwargs.get("minimize_phase_ramp",False)
    do_phase_match = kwargs.get("real_space_phase_match",True)
    do_align_com_support = kwargs.get("align_com_support",False)

    if logger0 != None:
        s = "  "
        if enantio: s += "Enantio, "
        if shifted: s += "Shifted, "
        if pixels_to_exclude != None: s += "%i pixels specified to exclude, " % pixel_to_exclude.sum()
        if do_maximize_overlap: s += "Maximizing overlap, "
        if do_minimize_phase_ramp: s += "Minimize phase ramp, "
        if do_phase_match: s += "Phase match, "
        if do_align_com_support: s += "Align center of mass of support, "
        if center_result != None: s += "Center result: %s, " % center_result
        logger0.debug("PRTF runs with the folloing configuration: %s",s[:-2])

    if debug:
        os.system("rm %s/testdata/prtf*" % this_folder)

    Nx = imgs0.shape[2]
    Ny = imgs0.shape[1]
    cx = kwargs.get("cx",(Nx-1)/2.)
    cy = kwargs.get("cy",(Ny-1)/2.)
    selection = kwargs.get("selection",numpy.ones(imgs0.shape[0],dtype="bool"))
    N = selection.sum()
    imgs = numpy.zeros(shape=(N,imgs0.shape[1],imgs0.shape[2]),dtype=imgs0.dtype)
    msks = numpy.zeros(shape=(N,msks0.shape[1],msks0.shape[2]),dtype="float")
    k = 0
    for i in range(imgs0.shape[0]):
        if selection[i]:
            if shifted:
                imgs[k,:,:] = numpy.fft.fftshift(imgs0[i,:,:])
                msks[k,:,:] = numpy.fft.fftshift(msks0[i,:,:])
            else:
                imgs[k,:,:] = imgs0[i,:,:]
                msks[k,:,:] = msks0[i,:,:]
            if debug and False:
                pylab.imsave(this_folder+"/testdata/prtf_%i_imgs%i.png" % (K,k),abs(imgs[k]))
                pylab.imsave(this_folder+"/testdata/prtf_%i_msks%i.png" % (K,k),abs(msks[k]))
                pylab.imsave(this_folder+"/testdata/prtf_%i_imgs0%i.png" % (K,k),abs(imgs0[i]))
                pylab.imsave(this_folder+"/testdata/prtf_%i_msks0%i.png" % (K,k),abs(msks0[i]))
            k += 1
    # Average reconstructions
    # superimpose for obtaining the averaged result of the reconstruction
    imgs1 = numpy.zeros(shape=(N,Ny,Nx),dtype=imgs0.dtype)
    msks1 = numpy.zeros(shape=(N,Ny,Nx),dtype="float")
    for i in range(0,N):
        img1 = imgs[i,:,:].copy()
        msk1 = msks[i,:,:].copy()
        img0 = imgs1[0,:,:].copy()
        msk0 = msks1[0,:,:].copy()

        if debug:
            j=0
            pylab.imsave(this_folder+"/testdata/prtf_%i_%i_I_%i.png" % (K,i,j),abs(img1),vmin=0.,vmax=3.)
            pylab.imsave(this_folder+"/testdata/prtf_%i_%i_IP_%i.png" % (K,i,j),numpy.angle(img1) % (2.*numpy.pi),vmin=0.,vmax=2.*numpy.pi)
            pylab.imsave(this_folder+"/testdata/prtf_%i_%i_M_%i.png" % (K,i,j),abs(msk1))
        if do_minimize_phase_ramp:
            [img1,translation] = minimize_phase_ramp(img1,shifted=False,periodic_boundary=True)
            msk1 = numpy.int16(abs(fourier_translation(msk1,translation)).round())
        
        
        if debug:
            j+=1
            pylab.imsave(this_folder+"/testdata/prtf_%i_%i_I_%i.png" % (K,i,j),abs(img1),vmin=0.,vmax=3.)
            pylab.imsave(this_folder+"/testdata/prtf_%i_%i_IP_%i.png" % (K,i,j),numpy.angle(img1) % (2.*numpy.pi),vmin=0.,vmax=2.*numpy.pi)
            pylab.imsave(this_folder+"/testdata/prtf_%i_%i_A_M_%i.png" % (K,i,j),abs(msk1))
        if do_maximize_overlap and i!=0:
            [img1,translation,turned] = maximize_overlap(img0,img1,enantio)
            if turned: msk1 = fft_turn180(msk1)
            msk1 = abs(fourier_translation(msk1,translation))

        if debug:
            j+=1
            pylab.imsave(this_folder+"/testdata/prtf_%i_%i_I_%i.png" % (K,i,j),abs(img1),vmin=0.,vmax=3.)
            pylab.imsave(this_folder+"/testdata/prtf_%i_%i_IP_%i.png" % (K,i,j),numpy.angle(img1) % (2.*numpy.pi),vmin=0.,vmax=2.*numpy.pi)
            pylab.imsave(this_folder+"/testdata/prtf_%i_%i_A_M_%i.png" % (K,i,j),abs(msk1))
        if do_phase_match and i!=0:
            weights = abs(img0)*abs(img1)
            img1 = abs(img1)*numpy.exp(1.j*(numpy.angle(img1)+phase_match(img0,img1,weights)))

        if debug:
            j+=1
            pylab.imsave(this_folder+"/testdata/prtf_%i_%i_I_%i.png" % (K,i,j),abs(img1),vmin=0.,vmax=3.)
            pylab.imsave(this_folder+"/testdata/prtf_%i_%i_IP_%i.png" % (K,i,j),numpy.angle(img1) % (2.*numpy.pi),vmin=0.,vmax=2.*numpy.pi)
            pylab.imsave(this_folder+"/testdata/prtf_%i_%i_A_M_%i.png" % (K,i,j),abs(msk1))
            print "Power: %f" % (abs(img1).sum()/(abs(img0).sum()+numpy.finfo("float32").eps))
            print "Avg. phase: %f,%f" % ((msk0*numpy.angle(img1)).mean(),(msk0*numpy.angle(img0)).mean())
            print "Diff: %f" % (abs(img1-img0).mean()/abs(img0).mean())
        imgs1[i,:,:] = img1[:,:]
        msks1[i,:,:] = numpy.int16(abs(msk1).round())[:,:]
    imgs1_super = imgs1.mean(0)
    msks1_super = msks1.mean(0)
    # Make PRTF
    # go to fourier space
    fimgs = numpy.zeros_like(imgs)
    fimgs1 = numpy.zeros_like(imgs)
    for i in range(N):
        fimgs[i,:,:] = numpy.fft.fftn(imgs[i,:,:])
        fimgs1[i,:,:] = numpy.fft.fftn(imgs1[i,:,:])
    # mask zeros
    PRTF = numpy.zeros_like(imgs)
    tmp = abs(fimgs1) != 0.
    if pixels_to_exclude != None:
        tmp *=  pixels_to_exclude
    PRTF[tmp] = fimgs1[tmp]/abs(fimgs1[tmp])
    PRTF = abs(PRTF.mean(0))
    PRTF[(fimgs == 0).sum(0) != 0] = 0.
    PRTF = numpy.array(PRTF,dtype="float32")

    if debug:
        pylab.imsave(this_folder+"/testdata/superI%i.png" % numpy.random.randint(1000),abs(imgs1_super),vmin=0,vmax=2.)
        pylab.imsave(this_folder+"/testdata/superM%i.png" % numpy.random.randint(1000),abs(msks1_super),vmin=0,vmax=1.)

    if center_result != None:
        if center_result == "image":
            CM = center_of_mass(abs(imgs1_super))
        elif center_result == "support_times_image":
            CM = center_of_mass(msks1_super*abs(imgs1_super))
        elif center_result == "support":
            CM = center_of_mass(msks1_super)
        imgs1_super = pixel_translation(imgs1_super,CM,3)
        msks1_super = pixel_translation(msks1_super,CM,3)

    if do_align_com_support:
        com = center_of_mass(msks1_super)
        if com[0] > 0 and com[1] > 0:
            imgs1_super = fft_turn180(imgs1_super)
            msks1_super = abs(fft_turn180(msks1_super))

    if shifted:
        imgs1_super = numpy.fft.fftshift(imgs1_super)
        msks1_super = numpy.fft.fftshift(msks1_super)
        for i in range(N):
            fimgs1[i,:,:] = numpy.fft.fftshift(fimgs1[i,:,:])
            imgs1[i,:,:] = numpy.fft.fftshift(imgs1[i,:,:])
            msks1[i,:,:] = numpy.fft.fftshift(msks1[i,:,:])
        PRTF = numpy.fft.fftshift(PRTF)

    msks1_super = numpy.int16(msks1_super)
    msks1 = numpy.int16(msks1)
    return [PRTF,imgs1_super,msks1_super,imgs1,msks1,fimgs1]

def half_period_resolution(PRTF,pixel_edge_length,detector_distance,wavelength,cx=None,cy=None):
    # angular average of PRTF
    [r,PRTFr] = radial_mean(PRTF,cx=cx,cy=cy,rout=True)
    dx = wavelength/2./(numpy.sin(numpy.arctan(r*pixel_edge_length/detector_distance))+numpy.finfo('float64').eps)
    success = PRTFr > (1./numpy.e)
    if success.sum() == len(success):
        i = -1
    elif success.sum() > 0:
        i = (numpy.arange(len(success))[success==False])[0]-1
    else:
        i = 0
    return [dx[i],PRTFr,1./dx]

# under testing
def half_pixel_downsampling(img0,msk0,downsampling=2,mode="conservative",cx=None,cy=None,cropLength=None):
    Nx = int(numpy.ceil(img0.shape[1]/(1.*downsampling)))*downsampling
    Ny = int(numpy.ceil(img0.shape[0]/(1.*downsampling)))*downsampling
    img = numpy.zeros(shape=(Ny,Nx),dtype="float")
    msk = numpy.zeros(shape=(Ny,Nx),dtype="int16")
    img[:img0.shape[0],:img0.shape[1]] = img0[:,:]
    msk[:msk0.shape[0],:msk0.shape[1]] = msk0[:,:]
    X,Y = numpy.meshgrid(numpy.arange(Nx),numpy.arange(Ny))
    Nx_XxX = Nx/downsampling + ((Nx%downsampling)!=0)
    Ny_XxX = Ny/downsampling + ((Ny%downsampling)!=0)
    X_XxX,Y_XxX = numpy.meshgrid(numpy.arange(Nx_XxX),numpy.arange(Ny_XxX))

    I_XxX = numpy.zeros(shape=(Ny,Nx),dtype="int")
    M = numpy.zeros(shape=(Ny,Nx),dtype="int")
    k = 0
    for iy_XxX in numpy.arange(Ny_XxX):
        for ix_XxX in numpy.arange(Nx_XxX):
            x0 = ix_XxX*downsampling
            y0 = iy_XxX*downsampling
            k += 1
            popbox = range(downsampling*downsampling)
            for i in numpy.arange(downsampling*downsampling):
                j = popbox.pop(numpy.random.randint(len(popbox)))
                x = j%downsampling
                y = j/downsampling
                if y0+y<Ny and x0+x<Nx:
                    M[y0+y,x0+x] = i
                    I_XxX[y0+y,x0+x] = k
                    

    ind0 = M<downsampling*downsampling/2
    ind1 = M>=downsampling*downsampling/2
    sind0 = I_XxX[ind0].argsort()#[::-1]
    sind1 = I_XxX[ind1].argsort()#[::-1]
    
    B0_img = numpy.reshape((img[ind0])[sind0],(Nx_XxX*Ny_XxX,downsampling*downsampling/2))
    B1_img = numpy.reshape((img[ind1])[sind1],(Nx_XxX*Ny_XxX,downsampling*downsampling/2))

    B0_msk = numpy.reshape((msk[ind0])[sind0],(Nx_XxX*Ny_XxX,downsampling*downsampling/2))
    B1_msk = numpy.reshape((msk[ind1])[sind1],(Nx_XxX*Ny_XxX,downsampling*downsampling/2))

    mskXxX0 = numpy.zeros(shape=(Ny_XxX,Nx_XxX),dtype="int16").flat
    mskXxX1 = numpy.zeros(shape=(Ny_XxX,Nx_XxX),dtype="int16").flat
    mskgoodXxX0 = numpy.zeros(shape=(Ny_XxX,Nx_XxX),dtype="int16").flat
    mskgoodXxX1 = numpy.zeros(shape=(Ny_XxX,Nx_XxX),dtype="int16").flat
    imgXxX0 = numpy.zeros(shape=(Ny_XxX,Nx_XxX),dtype="float").flat
    imgXxX1 = numpy.zeros(shape=(Ny_XxX,Nx_XxX),dtype="float").flat
    N0 = numpy.zeros(shape=(Ny_XxX,Nx_XxX),dtype="float").flat
    N1 = numpy.zeros(shape=(Ny_XxX,Nx_XxX),dtype="float").flat
    for i in range(downsampling*downsampling/2):
        m0 = B0_msk[:,i]
        m1 = B1_msk[:,i]
        n0 = ((m0 & cxitools.PIXEL_IS_IN_MASK) == 0)
        n1 = ((m1 & cxitools.PIXEL_IS_IN_MASK) == 0)
        i0 = B0_img[:,i]
        i1 = B1_img[:,i]
        mskXxX0 |= m0
        mskXxX1 |= m1
        mskgoodXxX0 |= n0*m0
        mskgoodXxX1 |= n1*m1
        imgXxX0 += n0*i0
        imgXxX1 += n1*i1
        N0 += n0
        N1 += n1
    if mode == "conservative":
        M = (N0+N1)==downsampling*downsampling
        imgXxX0[M==False] = 0
        imgXxX1[M==False] = 0
        imgXxX0[M] *= 2
        imgXxX1[M] *= 2
    elif mode == "non-conservative":
        M = (N0!=0)*(N1!=0)
        imgXxX0[M==False] = 0
        imgXxX1[M==False] = 0
        imgXxX0[M] *= downsampling*downsampling/N0[M]
        imgXxX1[M] *= downsampling*downsampling/N1[M]
        mskXxX0[M] = mskgoodXxX0[M]
        mskXxX1[M] = mskgoodXxX1[M]
    else:
        print "Invalid mode."
        return
        
    imgXxX0 = numpy.reshape(imgXxX0,(Ny_XxX,Nx_XxX))
    imgXxX1 = numpy.reshape(imgXxX1,(Ny_XxX,Nx_XxX))
    mskXxX0 = numpy.reshape(mskXxX0,(Ny_XxX,Nx_XxX))
    mskXxX1 = numpy.reshape(mskXxX1,(Ny_XxX,Nx_XxX))

    if cropLength != None:
        if cx == None or cy == None:
            center = "middle"
        else:
            center = [cy,cx]
        imgXxX0 = crop(imgXxX0,cropLength,center)
        imgXxX1 = crop(imgXxX1,cropLength,center)
        mskXxX0 = crop(mskXxX0,cropLength,center)
        mskXxX1 = crop(mskXxX1,cropLength,center)
        
    return [imgXxX0,mskXxX0,imgXxX1,mskXxX1]


def array_to_array(A1,A2,p0=None,origin="corner",mode="sum",fill_value=0.):
    N1 = numpy.array(A1.shape)
    N2 = numpy.array(A2.shape)
    d = len(N1)
    if d > 3:
        logger.error("Cannot handle more than 3 dimensional data.")
        return
    if p0 == None:
        p1 = numpy.zeros(d)
    else:
        p1 = p0
    if origin == "corner":
        p = p1
    elif origin == "middle":
        p = p1+(N2-1)/2.
    else:
        p = p1+origin
    p_min = numpy.int16((p-N1/2).round())
    p_max = p_min + N1
    print p_min,p_max
    A2_new = A2.copy()
    N2_new = N2.copy()
    origin_offset = numpy.zeros(d)
    for di in range(d):
        if p_min[di] < 0:
            offset = -p_min[di]
            N2_new[di] += offset
            p_min[di] += offset
            p_max[di] += offset
        if p_max[di] >= N2[di]:
            N2_new[di] = p_max[di]+1
    A2_new = numpy.zeros(shape=tuple(N2_new),dtype=A2.dtype) + fill_value
    print p_min,p_max
    if mode == "sum": f = lambda a,b: a+b
    elif mode == "replace": f = lambda a,b: b
    else: logger.error("%s is not a valid mode." % mode)
    if d == 1:
        A2_new[p_min[0]:p_max[0]] = f(A2_new[p_min[0]:p_max[0]],A1[:])
    elif d == 2:
        A2_new[p_min[0]:p_max[0],p_min[1]:p_max[1]] = f(A2_new[p_min[0]:p_max[0],p_min[1]:p_max[1]],A1[:,:])
    elif d == 3:
        A2_new[p_min[0]:p_max[0],p_min[1]:p_max[1],p_min[2]:p_max[2]] = f(A2_new[p_min[0]:p_max[0],p_min[1]:p_max[1],p_min[2]:p_max[2]],A1[:,:,:])
    return A2_new
        
def get_outline_pixels(M):
    import scipy.signal
    K = numpy.ones(shape=(3,3),dtype="int16")
    return (scipy.signal.convolve2d(M,K,mode='same',boundary='fill') != 0) * (M == 0)

def extend_mask_blurry(M,N):
    import scipy.signal
    Y,X = numpy.indices((N*2+1,N*2+1))
    Y -= N/2
    X -= N/2
    K = numpy.array((numpy.sqrt(X**2+Y**2)<=N),dtype="float")
    K /= K.sum()
    return scipy.signal.convolve2d(M,K,mode='same',boundary='fill')


