import logging
logger = logging.getLogger("Propagator")
import numpy
import config,proptools


def calculatePattern_SampleSphere(input_obj):
    # scattering amplitude from homogeneous sphere: F = sqrt(I_0 Omega_p) 2pi/wavelength^2 [ 4/3 pi R^3  3 { sin(qR) - qR cos(qR) } / (qR)^3 ] dn_real

    wavelength = input_obj.source.photon.get_wavelength()
    I_0 = input_obj.source.get_intensity("ph/m2")
    Omega_p = input_obj.detector.get_pixel_solid_angle("binned")
    if input_obj.sample.material == None:
        dn_real = 1.
    else:
        dn_real = input_obj.sample.material.get_dn().real

    R = input_obj.sample.radius
    V = 4/3.*numpy.pi*R**3
    K = I_0*Omega_p*(2*numpy.pi/wavelength**2*V*dn_real)**2

    q = input_obj.generate_absqmap()
    F = proptools.F_sphere_diffraction(K,q,R)
    return F

def calculatePattern_SampleSpheroid(input_obj):
    # scattering amplitude from homogeneous spheroid (small-angle approximation):
    # theta: rotation around x-axis (1st)
    # phi: rotation around (beam) z-axis (2nd)
    # a: semidiameter perpendicular to rotation axis
    # c: semidiameter parallel to rotation axis
    # F = sqrt(I_0 Omega_p) 2pi/wavelength^2 [ 4/3 pi a^2 c  3 { sin(qH) - qH cos(qH) } / (qH)^3 ] dn_real
    # H = sqrt(a^2 sin^2(g)+c^2 cos^2(g))
    # g = arccos( ( -qX sin(phi) cos(theta) + qY cos(phi) cos(theta) ) / sqrt(qX^2+qY^2) )

    wavelength = input_obj.source.photon.get_wavelength()
    I_0 = input_obj.source.pulse_energy / input_obj.source.photon.get_energy("J") / input_obj.source.get_area() # [I_0] = photons/m**2
    Omega_p = (input_obj.detector.pixel_size*input_obj.detector.binning)**2 / input_obj.detector.distance**2
    dn_real = (1-input_obj.sample.material.get_n()).real
    
    a = input_obj.sample.a
    c = input_obj.sample.c
    theta = input_obj.sample.theta
    phi = input_obj.sample.phi
    V = 4/3.*numpy.pi*a**2*c
    K = I_0*Omega_p*(2*numpy.pi/wavelength**2*V*dn_real)**2

    [qX,qY,qZ] = input_obj.generate_qmap() 

    F = proptools.F_spheroid_diffraction(K,qX,qY,a,c,theta,phi)

    return F

def calculatePattern_SampleSpheres(input_obj):
    wavelength = input_obj.source.photon.get_wavelength()
    I_0 = input_obj.source.pulse_energy / input_obj.source.photon._energy / input_obj.source.get_area() # [I_0] = photons/m**2
    Omega_p = (input_obj.detector.pixel_size*input_obj.detector.binning)**2 / input_obj.detector.distance**2
    radii = input_obj.sample.get_radii()
    dn_real = (1-input_obj.sample.material.get_n()).real
    absq = input_obj.generate_absqmap()
    q = input_obj.generate_qmap()
    F = numpy.zeros(shape=absq.shape,dtype='complex')
    for R in radii:
        Fr = numpy.zeros_like(F)
        Fr[absq!=0.0] = (numpy.sqrt(I_0*Omega_p)*2*numpy.pi/wavelength**2*4/3.0*numpy.pi*R**3*3*(numpy.sin(absq[absq!=0.0]*R)-absq[absq!=0.0]*R*numpy.cos(absq[absq!=0.0]*R))/(absq[absq!=0.0]*R)**3*dn_real)
        try: Fr[absq==0] = numpy.sqrt(I_0*Omega_p)*2*numpy.pi/wavelength**2*4/3.0*numpy.pi*R**3*dn_real
        except: pass

        indices = input_obj.sample.r==R

        for i in range(sum(indices)):
            looger.debug("%i" % i)
            d = [numpy.array(input_obj.sample.z)[indices][i],
                 numpy.array(input_obj.sample.y)[indices][i],
                 numpy.array(input_obj.sample.x)[indices][i]]
            F[:,:] += (Fr*input_obj.get_phase_ramp(q,d))[:,:]
    return F


def calculatePattern_SampleMap(input_obj):
    wavelength = input_obj.source.photon.get_wavelength()
    I_0 = input_obj.source.pulse_energy / input_obj.source.photon._energy / input_obj.source.get_area() # [I_0] = photons/m**2
    Omega_p = (input_obj.detector.get_pixel_size('binned'))**2 / input_obj.detector.distance**2
    # scattering amplitude from dn-map: F = sqrt(I_0 Omega_p) 2pi/wavelength^2 [ DFT{dn_perp} ] dA
    #print input_obj.sample.dX/input_obj.get_real_space_resolution_element()
    if abs(input_obj.sample.dX/input_obj.get_real_space_resolution_element()-1) < 0.05:
        pass
    elif input_obj.sample.dX <= input_obj.get_real_space_resolution_element():
        input_obj.sample.map3d_fine = input_obj.sample.map3d.copy()
        d = input_obj.get_real_space_resolution_element()/input_obj.sample.dX
        N_mapfine = input_obj.sample.map3d_fine.shape[0]
        N_map = int(round(N_mapfine/d))
        input_obj.sample.map3d = imgutils.lanczos_interp(input_obj.sample.map3d_fine,N_map)
        input_obj.sample.dX_fine = numpy.copy(input_obj.sample.dX)
        input_obj.sample.dX = N_mapfine/(1.*N_map)*input_obj.sample.dX
    else:
        if 'input_obj.sample.dX_fine' in locals():
            if input_obj.sample.dX_fine <= input_obj.get_real_space_resolution_element():
                d = input_obj.get_real_space_resolution_element()/input_obj.sample.dX
                N_mapfine = input_obj.sample.map3d_fine.shape[0]
                N_map = int(round(N_mapfine/d))
                input_obj.sample.map3d = imgutils.lanczos_interp(input_obj.sample.map3d_fine,N_map)
                input_obj.sample.dX = N_mapfine/(1.*N_map)*input_obj.sample.dX
            else:
                logger.error("Finer real space sampling required for chosen geometry.")
                return
        else:
            logger.error("Finer real space sampling required for chosen geometry.")
            return

    q_scaled = input_obj.generate_qmap(True)
    valid_mask = (abs(q_scaled[:,:,0])<0.5)*(abs(q_scaled[:,:,1])<0.5)*(abs(q_scaled[:,:,2])<0.5)
    qx = q_scaled[:,:,2]
    qy = q_scaled[:,:,1]
    qz = q_scaled[:,:,0]
    qx[valid_mask==False] = 0.
    qy[valid_mask==False] = 0.
    qz[valid_mask==False] = 0.

    logger.debug("Propagate pattern of %i x %i pixels." % (q_scaled.shape[1],q_scaled.shape[0]))
    F = numpy.sqrt(I_0*Omega_p)*2*numpy.pi/wavelength**2 \
        * xcorepropagation.nfftSingleCore(input_obj.sample.map3d,q_scaled) \
        * input_obj.sample.dX**3

    F *= valid_mask
    

    #F = numpy.sqrt(I_0*Omega_p)*2*numpy.pi/wavelength**2 \
    #    * xcorepropagation.nfftXCore(input_obj.sample.map3d-numpy.median(input_obj.sample.map3d),
    #                                 q_scaled,
    #                                 input_obj.propagation.N_processes) \
    #                                 * input_obj.sample.dX**3
    logger.debug("Got pattern of %i x %i pixels." % (F.shape[1],F.shape[0]))
    F = imgutils.crop(F,numpy.array([input_obj.detector.mask.shape[0],input_obj.detector.mask.shape[1]]))
    return F
