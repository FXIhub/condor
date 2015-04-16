#-----------------------------------------------------------------------------------------------------
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

import numpy, sys, numpy, types, pickle, time, math
import icosahedron,imgutils
 
import logging
logger = logging.getLogger("Condor")
from log import log

def random_euler_angles():
    """
    Generates a triplet (phi, theta, psi) of random Euler angles.
    """
    r1,r2,r3 = numpy.random.random(3)
    q1 = numpy.sqrt(1.0-r1)*numpy.sin(2.0*numpy.pi*r2)
    q2 = numpy.sqrt(1.0-r1)*numpy.cos(2.0*numpy.pi*r2)
    q3 = numpy.sqrt(r1)*numpy.sin(2.0*numpy.pi*r3)
    q4 = numpy.sqrt(r1)*numpy.cos(2.0*numpy.pi*r3)
    e1 = math.atan2(2.0*(q1*q2+q3*q4), 1.0-2.0*(q2**2+q3**2))
    e2 = math.asin(2.0*(q1*q3-q4*q2))
    e3 = math.atan2(2.0*(q1*q4+q2*q3), 1.0-2.0*(q3**2+q4**2))
    return (e1,e2,e3)

def get_max_crystallographic_resolution(wavelength,min_detector_center_edge_distance,detector_distance):
    """
    Returns crystallographic resolution (full-period resolution at the closest edge)
    """
    return wavelength/numpy.sin(numpy.arctan(min_detector_center_edge_distance/detector_distance))
    
def get_nyquist_pixel_size(detector_distance,wavelength,particle_area):
    """
    Returns size of one Nyquist pixel on the detector in meter.
    """
    particle_radius = numpy.sqrt(particle_area/numpy.pi)
    return detector_distance * wavelength / (2*particle_radius)

def print_material_xray_properties(wavelength,thickness=1.0E-06,**margs):
    #r_0 = constants.value("classical electron radius")
    h =  constants.h
    c =  constants.c
    qe = constants.e

    photon_energy_eV = h*c/wavelength/qe
    M = Material(photon_energy_eV,**margs)
    n = M.get_n(photon_energy_eV)
    print n
    f = M.get_f(photon_energy_eV)
    dn = n-1
    delta = -dn.real
    beta = -dn.imag
    phase_shift = 2*numpy.pi*thickness*delta/wavelength
    #mu_a = 2*re*wavelength*f2
    #s = 1/mu_a/n
    T = numpy.exp(-4*numpy.pi*beta/wavelength*thickness)
    print "BEAM:"
    print "Wavelength = %.2f nm ; Energy = %.0f eV" % (wavelength/1.0E-09,photon_energy_eV)
    print "SAMPLE DENSITY:"
    print "Mass denstity: %.3f mg/cm^3" % M.massdensity
    print "SAMPLE SCATTERING AND ABSORPTION PARAMETERS:"
    print "Scattering factor (real part) f1 = %e" % f.real
    print "Scattering factor (imag part) f2 = %e" % f.imag
    print "Refraction coefficient n = 1 - delta - i beta"
    print "delta = %f" % delta
    print "beta = %f" % beta
    print "Phaseshift / %.2f um = %f pi" % (thickness/1.0E-6,phase_shift/numpy.pi)
    #print "Atomic photoabsorption cross section: mu_a = %f re^2" % (mu_a/re**2)
    #print "Attenuation length (drop off to 1/e): s = %f um" % (s/1.0E-6)
    print "Transmission after %.2f um sample: T = %.1f percent " % (thickness/1.0E-06,T*100)
    #atomic photoabsorption cross section

def rotation(vector_or_matrix,E0,E1,E2):
    cos = numpy.cos
    sin = numpy.sin
    #Lsq = vector[0]**2+vector[1]**2+vector[2]**2
    M = numpy.array([[cos(E1)*cos(E2),
                      -cos(E0)*sin(E2)+sin(E0)*sin(E1)*cos(E2),
                      sin(E0)*sin(E2)+cos(E0)*sin(E1)*cos(E2)],
                     [cos(E1)*sin(E2),
                      cos(E0)*cos(E2)+sin(E0)*sin(E1)*sin(E2),
                      -sin(E0)*cos(E2)+cos(E0)*sin(E1)*sin(E2)],
                     [-sin(E1),
                      sin(E0)*cos(E1),
                      cos(E0)*cos(E1)]])
    rotated = numpy.dot(M,vector_or_matrix)
    #Lsq_rotated = rotated_vector[0]**2+rotated_vector[1]**2+rotated_vector[2]**2
    #if abs((Lsq - Lsq_rotated)/Lsq) > 0.001:
    #    print "ERROR: Rotation changes length!"
    #    print Lsq,Lsq_rotated
    return rotated

def generate_absqmap(X,Y,pixel_size,detector_distance,wavelength):
    qmap = generate_qmap(X,Y,pixel_size,detector_distance,wavelength)
    qmap = numpy.sqrt(qmap[:,:,0]**2+qmap[:,:,1]**2+qmap[:,:,2]**2)
    return qmap

def generate_qmap(X,Y,pixel_size,detector_distance,wavelength,euler_angle_0=0.,euler_angle_1=0.,euler_angle_2=0.):
    log(logger.debug,"Allocating qmap.")
    R_Ewald = 2*numpy.pi/wavelength
    qx = R_Ewald*2*numpy.sin(numpy.arctan2(pixel_size*X,detector_distance)/2.)
    qy = R_Ewald*2*numpy.sin(numpy.arctan2(pixel_size*Y,detector_distance)/2.)
    phi = numpy.arctan2(pixel_size*numpy.sqrt(X**2+Y**2),detector_distance)
    qz = R_Ewald*(1-numpy.cos(phi))
    qmap = numpy.zeros(shape=(X.shape[0],Y.shape[1],3))
    qmap[:,:,0] = qz[:,:]
    qmap[:,:,1] = qy[:,:]
    qmap[:,:,2] = qx[:,:]
    if euler_angle_0 != 0. or euler_angle_1 != 0. or euler_angle_2 != 0.:
        log(logger.debug,"Applying qmap rotation with angles %f, %f, %f." % (euler_angle_0,euler_angle_1,euler_angle_2))
        # Old and slow
        #for iy in numpy.arange(0,qmap.shape[0]):
        #    for ix in numpy.arange(0,qmap.shape[1]):
        #        qmap[iy,ix,:] = rotation(qmap[iy,ix,:],euler_angle_0,euler_angle_1,euler_angle_2)
        cos = numpy.cos
        sin = numpy.sin
        E0 = euler_angle_0
        E1 = euler_angle_1
        E2 = euler_angle_2
        M = numpy.array([[cos(E1)*cos(E2),
                          -cos(E0)*sin(E2)+sin(E0)*sin(E1)*cos(E2),
                          sin(E0)*sin(E2)+cos(E0)*sin(E1)*cos(E2)],
                         [cos(E1)*sin(E2),
                          cos(E0)*cos(E2)+sin(E0)*sin(E1)*sin(E2),
                          -sin(E0)*cos(E2)+cos(E0)*sin(E1)*sin(E2)],
                         [-sin(E1),
                          sin(E0)*cos(E1),
                          cos(E0)*cos(E1)]])
        Y,X = numpy.mgrid[:qmap.shape[0],:qmap.shape[1]]
        Y = Y.flatten()
        X = X.flatten()
        s = qmap.shape
        qmap = numpy.array([numpy.dot(M,qmap[iy,ix,:]) for ix,iy in zip(X,Y)])
        qmap = qmap.reshape(s)
    return qmap

def generate_qmap_ori(X,Y,pixel_size,detector_distance,wavelength):
    phi = numpy.arctan2(pixel_size*numpy.sqrt(X**2+Y**2),detector_distance)
    R_Ewald = 2*numpy.pi/wavelength
    qx = R_Ewald*2*numpy.sin(numpy.arctan2(pixel_size*X,detector_distance)/2.)
    qy = R_Ewald*2*numpy.sin(numpy.arctan2(pixel_size*Y,detector_distance)/2.)
    qz = R_Ewald*(1-numpy.cos(phi))
    qmap = numpy.zeros(shape=(X.shape[0],Y.shape[1],3))
    qmap[:,:,0] = qz[:,:]
    qmap[:,:,1] = qy[:,:]
    qmap[:,:,2] = qx[:,:]
    return qmap



x_to_q = lambda x,pixel_size,detector_distance,wavelength: 4*numpy.pi/wavelength*numpy.sin(numpy.arctan(x*pixel_size/detector_distance)/2.)

def get_sphere_diffraction_extrema_positions(N=1,xtol=1E-12,plot=False):
    from scipy.optimize import fmin

    f = lambda s:  numpy.sin(2*numpy.pi*s) - 2*numpy.pi*s*numpy.cos(2*numpy.pi*s)
    f_err = lambda s: f(s)**2
    g = lambda s:  2*numpy.pi*numpy.cos(2*numpy.pi*s) - 2*numpy.pi*numpy.cos(2*numpy.pi*s) + 4*numpy.pi**2*s*numpy.sin(2*numpy.pi*s)
    g_err = lambda s: g(s)**2

    s_mins = []
    s_maxs = []
    s_guess = 0.7
    while len(s_mins) < N:
        s_min = fmin(f_err,s_guess,(),xtol, ftol=0.0001, maxiter=None, maxfun=None, full_output=0, disp=0)
        s_max = fmin(g_err,s_min+0.2,(),xtol, ftol=0.0001, maxiter=None, maxfun=None, full_output=0, disp=0)
        s_guess = s_max+0.2
        s_mins.append(s_min)
        s_maxs.append(s_max)

    s_mins = numpy.array(s_mins)
    s_maxs = numpy.array(s_maxs)

    if plot:
        import matplotlib
        x = numpy.arange(0,s_mins.max()+1.,0.01)
        matplotlib.plot(x,f(x))
        matplotlib.plot(s_mins,numpy.zeros(len(s_mins)),".")
        matplotlib.plot(s_maxs,numpy.zeros(len(s_maxs)),".")
        matplotlib.show()

    return [s_mins,s_maxs]

def mask_fringe_sphere(q,r,i_fringe):
    s_mins = get_sphere_diffraction_extrema_positions(i_fringe+1)[0]
    s = q*r/2./numpy.pi
    M = s < s_mins[i_fringe]
    if i_fringe != 0:
        M *= s > s_mins[i_fringe-1]
    return M

def mask_min_sphere(q,r,i_min,ds=0.25):
    s_mins = get_sphere_diffraction_extrema_positions(i_min+1)[0]
    s = q*r/2./numpy.pi
    M = (s > s_mins[i_min]-ds)*(s < s_mins[i_min]+ds)
    return M

def mask_fringe_spheroid(qX,qY,a,c,theta,phi,i_fringe):
    s_mins = get_sphere_diffraction_extrema_positions(i_fringe+1)[0]
    H = _H_spheroid_diffraction(qX,qY,a,c,theta,phi)
    q = numpy.sqrt(qX**2+qY**2)
    s = q*H/2./numpy.pi
    M = s < s_mins[i_fringe]
    if i_fringe != 0:
        M *= s > s_mins[i_fringe-1]
    return M

def mask_min_spheroid(qX,qY,a,c,theta,phi,i_min,ds=0.25):
    s_mins = get_sphere_diffraction_extrema_positions(i_min+1)[0]
    H = _H_spheroid_diffraction(qX,qY,a,c,theta,phi)
    q = numpy.sqrt(qX**2+qY**2)
    s = q*H/2./numpy.pi
    M = (s > s_mins[i_min]-ds)*(s < s_mins[i_min]+ds)
    return M

def spheroid_support(N,dx,a,c,phi):
    Y,X = numpy.indices((N,N))
    X = numpy.float64(X-N/2)*dx
    Y = numpy.float64(Y-N/2)*dx
    T = numpy.arctan2(Y,X)-phi
    R = numpy.sqrt(X**2+Y**2)
    Rmax = a*c/numpy.sqrt((a*numpy.sin(T))**2+(c*numpy.cos(T))**2)
    M = R<=Rmax
    return M
                


# scattering amplitude from homogeneous sphere:
# -----------------------------------------------
# Source:
# Feigin 1987
#
# r: sphere radius
#
# F = sqrt(I_0) rho_e p/D r_0 4/3 pi r^3 [ 3 { sin(qr) - qr cos(qr) } / (qr)^3 ]
#   = sqrt(I_0) rho_e p/D r_0 V f(r,qx,qy)
# f = 3 { sin(qr) - qr cos(qr) } / (qr)^3
# K = I_0 (rho_e p/D r_0 V)^2
# S = I_0 rho_e^2 = K / (p/D r_0 V)^2
# ============================================================================================
# I = F^2 = K [ f(r,qx,qy) ]^2
# ============================================================================================
_F_sphere_diffraction = lambda K,q,r: numpy.sqrt(abs(K))*3*(numpy.sin(q*r)-q*r*numpy.cos(q*r))/((q*r)**3+numpy.finfo("float64").eps)
F_sphere_diffraction = lambda K,q,r: ((q*r)**6 < numpy.finfo("float64").resolution)*numpy.sqrt(abs(K)) + ((q*r)**6 >= numpy.finfo("float64").resolution)*_F_sphere_diffraction(K,q,r)
_I_sphere_diffraction = lambda K,q,r: abs(K)*(3*(numpy.sin(q*r)-q*r*numpy.cos(q*r))/((q*r)**3+numpy.finfo("float64").eps))**2
I_sphere_diffraction = lambda K,q,r: ((q*r)**6 < numpy.finfo("float64").resolution)*abs(K) + ((q*r)**6 >= numpy.finfo("float64").resolution)*_I_sphere_diffraction(K,q,r)
Fringe_sphere_diffraction = None

# scattering amplitude from homogeneous spheroid:
# -----------------------------------------------
# Sources:
# Feigin 1987
# Hamzeh,Bragg 1974
#
# a: radius perpendicular to the rotation axis of the ellipsoid (a)
# c: radius along the rotation axis of the ellipsoid (c)
# theta: rotation around x-axis (1st)
# phi: rotation around z-axis (2nd)
#
# F = sqrt(I_0) rho_e p/D r_0 4/3 pi a^2 c [ 3 { sin(qH) - qH cos(qH) } / (qH)^3 ]
#   = sqrt(I_0) rho_e p/D r_0 V f(a,c,theta,phi,qx,qy)
# f = 3 { sin(qH) - qH cos(qH) } / (qH)^3
# H = sqrt(asq sin^2(g)+csq cos^2(g))
# g = arccos( ( -qX cos(theta) sin(phi) + qY cos(theta) cos(phi) ) / sqrt(qX^2+qY^2) )
# K = I_0 (rho_e p/D r_0 V)^2
# S = I_0 rho_e^2 = K / (p/D r_0 V)^2
# ============================================================================================
# I = K [ f(asq,csq,theta,phi,qx,qy) ]^2
# ============================================================================================
_q_spheroid_diffraction = lambda qX,qY: numpy.sqrt(qX**2+qY**2)
_g_spheroid_diffraction = lambda qX,qY,theta,phi: numpy.arccos((-qX*numpy.cos(theta)*numpy.sin(phi)+qY*numpy.cos(theta)*numpy.cos(phi))/(_q_spheroid_diffraction(qX,qY)+numpy.finfo("float64").eps))
_H_spheroid_diffraction = lambda qX,qY,a,c,theta,phi: numpy.sqrt(a**2*numpy.sin(_g_spheroid_diffraction(qX,qY,theta,phi))**2+c**2*numpy.cos(_g_spheroid_diffraction(qX,qY,theta,phi))**2)
_qH_spheroid_diffraction = lambda qX,qY,a,c,theta,phi: _q_spheroid_diffraction(qX,qY)*_H_spheroid_diffraction(qX,qY,a,c,theta,phi)
_F_spheroid_diffraction = lambda K,qX,qY,a,c,theta,phi: numpy.sqrt(abs(K))*3.*(numpy.sin(_qH_spheroid_diffraction(qX,qY,a,c,theta,phi))-_qH_spheroid_diffraction(qX,qY,a,c,theta,phi)*numpy.cos(_qH_spheroid_diffraction(qX,qY,a,c,theta,phi)))/(_qH_spheroid_diffraction(qX,qY,a,c,theta,phi)**3+numpy.finfo("float64").eps)
F_spheroid_diffraction = lambda K,qX,qY,a,c,theta,phi: (_qH_spheroid_diffraction(qX,qY,a,c,theta,phi)**6 < numpy.finfo("float64").resolution)*numpy.sqrt(abs(K)) + (_qH_spheroid_diffraction(qX,qY,a,c,theta,phi)**6 >= numpy.finfo("float64").resolution)*_F_spheroid_diffraction(K,qX,qY,a,c,theta,phi)
_I_spheroid_diffraction = lambda K,qX,qY,a,c,theta,phi: abs(K)*(3.*(numpy.sin(_qH_spheroid_diffraction(qX,qY,a,c,theta,phi))-_qH_spheroid_diffraction(qX,qY,a,c,theta,phi)*numpy.cos(_qH_spheroid_diffraction(qX,qY,a,c,theta,phi)))/(_qH_spheroid_diffraction(qX,qY,a,c,theta,phi)**3+numpy.finfo("float64").eps))**2
I_spheroid_diffraction = lambda K,qX,qY,a,c,theta,phi: (_qH_spheroid_diffraction(qX,qY,a,c,theta,phi)**6 < numpy.finfo("float64").resolution)*abs(K) + (_qH_spheroid_diffraction(qX,qY,a,c,theta,phi)**6 >= numpy.finfo("float64").resolution)*_I_spheroid_diffraction(K,qX,qY,a,c,theta,phi)

def get_spheroid_diffraction_formula(p,D,wavelength,X=None,Y=None):
    if X != None and Y != None:
        qmap = generate_qmap(X,Y,p,D,wavelength)
        qX = qmap[:,:,2]
        qY = qmap[:,:,1]
        qZ = qmap[:,:,0]
        I = lambda K,a,c,theta,phi: I_spheroid_diffraction(K,qX,qY,a,c,theta,phi)
    else:
        qmap = lambda X,Y: generate_qmap(X,Y,p,D,wavelength)
        I = lambda X,Y,K,a,c,theta,phi: I_sphere_diffraction(K,qmap(X,Y)[0],qmap(X,Y)[1],a,c,theta,phi)
    return I

def get_sphere_diffraction_formula(p,D,wavelength,X=None,Y=None):
    if X != None and Y != None:
        q = generate_absqmap(X,Y,p,D,wavelength)
        I = lambda K,r: I_sphere_diffraction(K,q,r)
    else:
        q = lambda X,Y: generate_absqmap(X,Y,p,D,wavelength)
        I = lambda X,Y,K,r: I_sphere_diffraction(K,q(X,Y),r)
    return I



to_spheroid_semi_diameter_a = lambda diameter,flattening: flattening**(1/3.)*diameter/2.
to_spheroid_semi_diameter_c = lambda diameter,flattening: flattening**(-2/3.)*diameter/2.
to_spheroid_diameter = lambda a,c: 2*(a**2*c)**(1/3.)
to_spheroid_flattening = lambda a,c: a/c

def to_spheroid_theta(euler_angle_0,euler_angle_1,euler_angle_2):
    v_z = numpy.array([1.0,0.0,0.0])
    v_y = numpy.array([0.0,1.0,0.0])
    v_rot = rotation(v_y,euler_angle_0,euler_angle_1,euler_angle_2)
    theta = numpy.arcsin(numpy.dot(v_rot,v_z))
    return theta

def to_spheroid_phi(euler_angle_0,euler_angle_1,euler_angle_2):
    v_y = numpy.array([0.0,1.0,0.0])
    v_rot = rotation(v_y,euler_angle_0,euler_angle_1,euler_angle_2)
    v_rot[0] = 0.0
    v_rot = v_rot / numpy.sqrt(v_rot[0]**2+v_rot[1]**2+v_rot[2]**2)       
    phi = numpy.arccos(numpy.dot(v_rot,v_y))
    return phi



def make_icosahedron_map(N,nRmax,euler1=0.,euler2=0.,euler3=0.):
    logger.debug("Building icosahedral geometry")
    logger.debug("Grid: %i x %i x %i (%i voxels)" % (N,N,N,N**3))
    t0 = time.time()
    icomap = icosahedron.icosahedron(N,nRmax,(euler1,euler2,euler3))
    t1 = time.time()
    logger.debug("Built map within %f seconds." % (t1-t0))
    return icomap

def make_icosahedron_map_old(N,nRmax,euler1=0.,euler2=0.,euler3=0.):
    na = nRmax/numpy.sqrt(10.0+2*numpy.sqrt(5))*4.
    nRmin = numpy.sqrt(3)/12*(3.0+numpy.sqrt(5))*na # radius at faces
    logger.debug("Building icosahedral geometry")
    n_list = imgutils.get_icosahedron_normal_vectors(euler1,euler2,euler3)
    X,Y,Z = 1.0*numpy.mgrid[0:N,0:N,0:N]
    X = X - (N-1)/2.
    Y = Y - (N-1)/2.
    Z = Z - (N-1)/2.
    logger.debug("Grid: %i x %i x %i (%i voxels)" % (N,N,N,N**3))
    icomap = numpy.zeros((len(n_list),N,N,N))
    # calculate distance of all voxels to all faces (negative inside, positive outside icosahedron)
    for i in range(len(n_list)):
        icomap[i,:,:,:] = (X*n_list[i][2]+Y*n_list[i][1]+Z*n_list[i][0])+nRmin
    s = 1.
    M = icomap.copy()
    temp = abs(M)<0.5*s
    icomap[temp] = 0.5+icomap[temp]/s
    icomap[M<(-0.5)*s] = 0
    icomap[M>0.5*s] = 1
    icomap = icomap.min(0)
    return icomap

def make_spheroid_map(N,nA,nB,euler0=0.,euler1=0.,euler2=0.):
    X,Y,Z = 1.0*numpy.mgrid[0:N,0:N,0:N]
    X = X-(N-1)/2.
    Y = Y-(N-1)/2.
    Z = Z-(N-1)/2.
    R_sq = X**2+Y**2+Z**2
    e_c = rotation(numpy.array([0.0,1.0,0.0]),euler0,euler1,euler2)
    d_sq_c = ((Z*e_c[0])+(Y*e_c[1])+(X*e_c[2]))**2
    r_sq_c = abs( R_sq * (1 - (d_sq_c/(R_sq+numpy.finfo("float32").eps))))
    spheroidmap = r_sq_c/nA**2+d_sq_c/nB**2
    spheroidmap[spheroidmap<=1] = 1
    spheroidmap[spheroidmap>1] = 0
    return spheroidmap

def make_sphere_map(N,nR):
    X,Y,Z = 1.0*numpy.mgrid[0:N,0:N,0:N]
    X = X-(N-1)/2.
    Y = Y-(N-1)/2.
    Z = Z-(N-1)/2.
    R = numpy.sqrt(X**2+Y**2+Z**2)
    spheremap = numpy.zeros(shape=R.shape,dtype="float64")
    spheremap[R<=nR] = 1
    spheremap[abs(nR-R)<0.5] = 0.5+0.5*(nR-R[abs(nR-R)<0.5])
    return spheremap


# Position conversion for a downsampled / upsampled array:
# downsample_pos
# pos: position in current binsize units 
# size: current array size
# binning: downsampled binsize / current binsize
downsample_pos = lambda pos,size,binning: (pos-(binning-1)/2.)*(size/(1.*binning)-1)/(1.*(size-binning))
# upsample_pos
# pos: position in current binsize units
# size: current array size
# binning: current binsize / upsampled binsize
upsample_pos   = lambda pos,size,binning: pos*(size*binning-binning)/(1.*(size-1))+(binning-1)/2.

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
        Nx_new = int(numpy.ceil(Nx/float(factor)))
        Ny_new = int(numpy.ceil(Ny/float(factor)))
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



def check_input(keys,req_keys,opt_keys,verbose=False):
    missing_keys = [k for k in req_keys if not isinstance(k,list) and (k not in keys)]
    for l in [l for k in keys if isinstance(k,list)]:
        tmp_miss = []
        # Alternatives (non-exclusive)
        for a in l:
            if isinstance(a,list):
                # Combined requirement
                com = [c for c in a if c not in keys]
                if len(com) > 0:
                    tmp_miss.append(com)
            else:
                # Single requirement
                if a not in keys:
                    tmp_miss.append(a)
        if len(tmp_miss) == len(l):
            missing_keys.append(tmp_miss)
    all_keys = req_keys + opt_keys
    def sublist_elements(l):
        l_new = []
        for k0 in l:
            if isinstance(k0,list):
                for k1 in k0:
                    l_new.append(k1)
            else:
                l_new.append(k0)
        return l_new
    illegal_keys = [k for k in keys if k not in all_keys]
    illegal_keys = [k for k in illegal_keys if k not in sublist_elements(all_keys)]
    illegal_keys = [k for k in illegal_keys if k not in sublist_elements(sublist_elements(all_keys))]
    if verbose:
        for illegal_key in illegal_keys:
            print "Illegal key: %s" % illegal_key
        if len(missing_keys) > 0:
            print "Missing key(s):"
        for missing_key in missing_keys:
            if isinstance(missing_key,list):
                print "= Alternatives:"
                for missing_key_alternative in missing_key:
                    if isinstance(missing_key_alternative,list):
                        s = "- ["
                        for missing_key_alternative_component in missing_key_alternative:
                            s += missing_key_alternative_component + " + "
                        s += "]"
                        print s
                    else:
                        print ("- " + missing_key_alternative)
            else:
                print ("= " + missing_key)
    return missing_keys,illegal_keys


gaussian = lambda x, sigma: numpy.exp(-x**2/(2*sigma**2))

gaussian_2dnorm = lambda x, sigma: gaussian(x, sigma) / ( 2 * numpy.pi * sigma**2 )


lorentzian = lambda x, sigma: sigma**2 / (x**2 + sigma**2)


_pseudo_lorenzian_A1 = 0.74447313315648778 
_pseudo_lorenzian_A2 = 0.22788162774723308
_pseudo_lorenzian_s1 = 0.73985516665883544
_pseudo_lorenzian_s2 = 2.5588165723260907
pseudo_lorentzian = lambda x, sigma: _pseudo_lorenzian_A1 * gaussian(x, _pseudo_lorenzian_s1*sigma) + \
                                     _pseudo_lorenzian_A2 * gaussian(x, _pseudo_lorenzian_s2*sigma)

pseudo_lorentzian_2dnorm = lambda x, sigma: pseudo_lorentzian(x, sigma) / ( 2. * numpy.pi * ( _pseudo_lorenzian_A1 * (_pseudo_lorenzian_s1*sigma)**2 + \
                                                                                              _pseudo_lorenzian_A2 * (_pseudo_lorenzian_s2*sigma)**2 ) )
