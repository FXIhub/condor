import numpy, matplotlib, sys, numpy, types, pickle, time, math

def random_euler_angles():
    """
    Generates a triplet (phi, theta, psi) of random Euler angles.
    """
    r1,r2,r3 = numpy.random(3)
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
    re = DICT_physical_constants['re']
    h = phy.DICT_physical_constants['h']
    c = phy.DICT_physical_constants['c']
    qe = phy.DICT_physical_constants['e']

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

def rotation(vector,phi,theta,psi):
    cos = numpy.cos
    sin = numpy.sin
    #Lsq = vector[0]**2+vector[1]**2+vector[2]**2
    M = numpy.array([[cos(theta)*cos(psi),-cos(phi)*sin(psi)+sin(phi)*sin(theta)*cos(psi),sin(phi)*sin(psi)+cos(phi)*sin(theta)*cos(psi)],
                     [cos(theta)*sin(psi),cos(phi)*cos(psi)+sin(phi)*sin(theta)*sin(psi),-sin(phi)*cos(psi)+cos(phi)*sin(theta)*sin(psi)],
                     [-sin(theta),sin(phi)*cos(theta),cos(phi)*cos(theta)]])
    rotated_vector = numpy.dot(M,vector)
    #Lsq_rotated = rotated_vector[0]**2+rotated_vector[1]**2+rotated_vector[2]**2
    #if abs((Lsq - Lsq_rotated)/Lsq) > 0.001:
    #    print "ERROR: Rotation changes length!"
    #    print Lsq,Lsq_rotated
    return rotated_vector

def generate_absqmap(X,Y,pixel_size,detector_distance,wavelength):
    qmap = generate_qmap(X,Y,pixel_size,detector_distance,wavelength)
    qmap = numpy.sqrt(qmap[:,:,0]**2+qmap[:,:,1]**2+qmap[:,:,2]**2)
    return qmap

def generate_qmap(X,Y,pixel_size,detector_distance,wavelength,euler_angle_0=0.,euler_angle_1=0.,euler_angle_2=0.):
    phi = numpy.arctan2(pixel_size*numpy.sqrt(X**2+Y**2),detector_distance)
    R_Ewald = 2*numpy.pi/wavelength
    qx = R_Ewald*2*numpy.sin(numpy.arctan2(pixel_size*X,detector_distance)/2.)
    qy = R_Ewald*2*numpy.sin(numpy.arctan2(pixel_size*Y,detector_distance)/2.)
    qz = R_Ewald*(1-numpy.cos(phi))
    qmap = numpy.zeros(shape=(X.shape[0],Y.shape[1],3))
    qmap[:,:,0] = qz[:,:]
    qmap[:,:,1] = qy[:,:]
    qmap[:,:,2] = qx[:,:]
    if euler_angle_0 != 0. or euler_angle_1 != 0. or euler_angle_2 != 0.:
        E0 = euler_angle_0
        E1 = euler_angle_1
        E2 = euler_angle_2
        M = numpy.array([[numpy.cos(E1)*numpy.cos(E2),
                          -numpy.cos(E0)*numpy.sin(E2)+numpy.sin(E0)*numpy.sin(E1)*numpy.cos(E2),
                          numpy.sin(E0)*numpy.sin(E2)+numpy.cos(E0)*numpy.sin(E1)*numpy.cos(E2)],
                         [numpy.cos(E1)*numpy.sin(E2),
                          numpy.cos(E0)*numpy.cos(E2)+numpy.sin(E0)*numpy.sin(E1)*numpy.sin(E2),
                          -numpy.sin(E0)*numpy.cos(E2)+numpy.cos(E0)*numpy.sin(E1)*numpy.sin(E2)],
                         [-numpy.sin(E1),
                          numpy.sin(E0)*numpy.cos(E1),
                          numpy.cos(E0)*numpy.cos(E1)]])
        for iy in numpy.arange(0,qmap.shape[0]):
            for ix in numpy.arange(0,qmap.shape[1]):
                qmap[iy,ix,:] = numpy.dot(M,qmap[iy,ix,:])
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
