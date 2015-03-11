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

import numpy

import logging
logger = logging.getLogger("Condor")

import utils
import condortools
from sample_species import SampleSpeciesSphere,SampleSpeciesSpheroid,SampleSpeciesMap

class Propagator:
    def __init__(self,source,sample,detector):
        self.source   = source
        self.sample   = sample
        self.detector = detector
        self._qmap_cache = {}
        self._map3d_cache = {}

    def propagate(self,**kwargs):
        N = self.sample.number_of_images
        O = {"source":{},"sample":{},"detector":{}}
        O_particles = {}
        N_particles_max = 0

        for i in range(N):
            logger.info("Calculation diffraction pattern (%i/%i). (PROGSTAT)" % (i+1,N))
            Os = self._propagate_single(**kwargs)
            for k in [k for k in Os.keys() if k not in ["source","sample","detector"]]:
                if k not in O:
                    O[k] = [None] * N
                O[k][i] = Os[k]
            # Source
            for k,v in Os["source"].items():
                if k not in O["source"]:
                    O["source"][k] = [None] * N
                O["source"][k][i] = Os["source"][k]
            # Sample
            for k in [k for k in Os["sample"].keys() if k != "particles"]:
                if k not in O["sample"]:
                    O["sample"][k] = [None] * N
                O["sample"][k][i] = Os["sample"][k]
            for j,p in zip(range(len(Os["sample"]["particles"])),Os["sample"]["particles"]):
                for k,v in p.items():
                    if k not in O_particles:
                        O_particles[k] = [None] * N
                    if O_particles[k][i] is None:
                        O_particles[k][i] = []
                    O_particles[k][i].append(v)
            N_particles_max = max([N_particles_max,len(Os["sample"]["particles"])])
            # Detector
            for k,v in Os["detector"].items():
                if k not in O["detector"]:
                    O["detector"][k] = [None] * N
                O["detector"][k][i] = Os["detector"][k]
        # Make arrays for particle parameters
        for k,v in O_particles.items():
            s = [N,N_particles_max]
            d = len(list(numpy.array(v).shape))
            for id in range(2,d):
                s.append(numpy.array(v).shape[id])
            s = tuple(s)
            A = numpy.zeros(shape=s,dtype=numpy.array(v[0]).dtype)
            for i in range(len(v)):

                A[i,:len(v[i])] = numpy.array(v[i])[:]
            O["sample"][k] = A 
        for k in O["source"].keys(): O["source"][k] = numpy.array(O["source"][k])
        for k in O["detector"].keys(): O["detector"][k] = numpy.array(O["detector"][k])
        for k in [k for k in O.keys() if k not in ["source","sample","detector"]]:
            O[k] = numpy.array(O[k])

        return O

    def _propagate_single(self,**kwargs):
        save_map = kwargs.get("save_map",False)
        save_qmap = kwargs.get("save_qmap",False)

        # Iterate objects
        D_source   = self.source.get_next()
        D_sample   = self.sample.get_next()
        D_detector = self.detector.get_next()

        # Pull out variables
        nx                  = D_detector["nx_binned"]
        ny                  = D_detector["ny_binned"]
        cx                  = D_detector["cx_binned"]
        cy                  = D_detector["cy_binned"]
        cx_unbinned         = D_detector["cx_unbinned"]
        cy_unbinned         = D_detector["cy_unbinned"]
        pixel_size          = D_detector["pixel_size_binned"]
        pixel_size_unbinned = D_detector["pixel_size_unbinned"]
        detector_distance   = D_detector["distance"]
        Omega_p             = D_detector["solid_angle_binned_pixel"]
        wavelength          = D_source["wavelength"]
        I_0                 = D_source["intensity_ph_m2"]

        # Calculate primary wave amplitude
        # F0 = sqrt(I_0 Omega_p) 2pi/wavelength^2
        F0 = numpy.sqrt(I_0*Omega_p)*2*numpy.pi/wavelength**2
        D_sample["F0"] = F0

        F_singles = []
        # Calculate patterns of all single particles individually
        for D_particle in D_sample["particles"]:
            p  = D_particle["_class_instance"]
            dn = p.material.get_dn(wavelength)
            e0 = D_particle["euler_angle_0"]
            e1 = D_particle["euler_angle_1"]
            e2 = D_particle["euler_angle_2"]
            if isinstance(p,SampleSpeciesSphere):
                # Scattering vectors
                qmap = self.get_qmap(nx=nx, ny=ny, cx=cx, cy=cy, pixel_size=pixel_size, detector_distance=detector_distance, wavelength=wavelength, 
                                     euler_angle_0=0., euler_angle_1=0., euler_angle_2=0.)
                q = numpy.sqrt(qmap[:,:,1]**2+qmap[:,:,2]**2)
                # Intensity scaling factor
                R = D_particle["diameter"]/2.
                V = 4/3.*numpy.pi*R**3
                K = (F0*V*dn.real)**2
                # Pattern
                F = condortools.F_sphere_diffraction(K,q,R)
            if isinstance(p,SampleSpeciesSpheroid):
                # Scattering vectors
                qmap = self.get_qmap(nx=nx, ny=ny, cx=cx, cy=cy, pixel_size=pixel_size, detector_distance=detector_distance, wavelength=wavelength, 
                                     euler_angle_0=0., euler_angle_1=0., euler_angle_2=0.)
                qx = qmap[:,:,2]
                qy = qmap[:,:,1]
                # Intensity scaling factor
                R = D_particle["diameter"]/2.
                V = 4/3.*numpy.pi*R**3
                K = (F0*V*dn.real)**2
                a = condortools.to_spheroid_semi_diameter_a(D_particle["diameter"],D_particle["flattening"])
                c = condortools.to_spheroid_semi_diameter_c(D_particle["diameter"],D_particle["flattening"])
                # Pattern
                # (Rotation is taken into account directly in the diffraction formula (much faster))
                theta = condortools.to_spheroid_theta(euler_angle_0=e0,euler_angle_1=e1,euler_angle_2=e2)
                phi = condortools.to_spheroid_phi(euler_angle_0=e0,euler_angle_1=e1,euler_angle_2=e2)
                F = condortools.F_spheroid_diffraction(K,qx,qy,a,c,theta,phi)
            if isinstance(p,SampleSpeciesMap):
                # Scattering vectors
                qmap = self.get_qmap(nx=nx, ny=ny, cx=cx, cy=cy, pixel_size=pixel_size, detector_distance=detector_distance, wavelength=wavelength,
                                     euler_angle_0=e0, euler_angle_1=e1, euler_angle_2=e2)
                # Intensity scaling factor
                R = D_particle["diameter"]/2.
                V = 4/3.*numpy.pi*R**3
                K = (F0*V*dn.real)**2
                # Generate map
                dx_required  = self.detector.get_real_space_resolution_element(wavelength,cx_unbinned,cy_unbinned) / numpy.sqrt(2)
                dx_suggested = self.detector.get_real_space_resolution_element_min(wavelength,cx_unbinned,cy_unbinned) / numpy.sqrt(2)
                map3d, dx = p.get_map3d(D_particle,dx_required,dx_suggested)
                dn_map3d = numpy.array(map3d,dtype="complex128") * dn
                if save_map:
                    D_particle["map3d"] = map3d
                    D_particle["dx"] = dx
                    D_particle["dx3"] = dx**3
                # Rescale and shape qmap for nfft
                qmap_scaled = dx * qmap / (2 * numpy.pi)
                #qmap /= self.get_absq_max(wavelength)/0.5*numpy.sqrt(2)
                qmap_shaped = qmap_scaled.reshape(qmap_scaled.shape[0]*qmap_scaled.shape[1],3)
                if save_qmap:
                    D_particle["qmap"] = qmap_shaped
                    #D_particle["qmap3d"] = detector.generate_qmap_ori(nfft_scaled=True)
                # Check inputs
                invalid_mask = (abs(qmap_shaped)>0.5)
                if (invalid_mask).sum() > 0:
                    qmap_shaped[invalid_mask] = 0.
                    logger.debug("%i invalid pixel positions." % invalid_mask.sum())
                logger.debug("Map3d input shape: (%i,%i,%i), number of dimensions: %i, sum %f" % (dn_map3d.shape[0],dn_map3d.shape[1],dn_map3d.shape[2],len(list(dn_map3d.shape)),abs(dn_map3d).sum()))
                if (numpy.isfinite(dn_map3d)==False).sum() > 0:
                    logger.warning("There are infinite values in the map3d of the object.")
                logger.debug("Scattering vectors shape: (%i,%i); Number of dimensions: %i" % (qmap_shaped.shape[0],qmap_shaped.shape[1],len(list(qmap_shaped.shape))))
                if (numpy.isfinite(qmap_shaped)==False).sum() > 0:
                    logger.warning("There are infinite values in the scattering vectors.")
                # NFFT
                fourier_pattern = utils.nfft.nfft(dn_map3d,qmap_shaped)
                # Check output - masking in case of invalid values
                if (invalid_mask).sum() > 0:
                    fourier_pattern[numpy.any(invalid_mask)] = numpy.nan
                # reshaping
                fourier_pattern = numpy.reshape(fourier_pattern,(qmap_scaled.shape[0],qmap_scaled.shape[1]))
                logger.debug("Got pattern of %i x %i pixels." % (fourier_pattern.shape[1],fourier_pattern.shape[0]))
                F = F0 * fourier_pattern * dx**3
            F_singles.append(F)

        F_tot = numpy.zeros_like(F)
        # Superimpose patterns
        for D_particle,F in zip(D_sample["particles"],F_singles):
            v = D_particle["position"]
            F_tot = F_tot + F * numpy.exp(-1.j*(v[0]*qmap[:,:,0]+v[1]*qmap[:,:,1]+v[2]*qmap[:,:,2])) 
        I_tot = self.detector.detect_photons(abs(F_tot)**2)
        M_tot = self.detector.get_mask(I_tot,False)
        M_tot_bitmask = self.detector.get_mask(I_tot,True)

        return {"fourier_pattern":F_tot,"intensity_pattern":I_tot,"mask_binary":M_tot,"mask":M_tot_bitmask,"source":D_source,"sample":D_sample,"detector":D_detector}
        
    def get_qmap(self, nx, ny, cx, cy, pixel_size, detector_distance, wavelength, euler_angle_0, euler_angle_1, euler_angle_2):
        calculate = False
        keys = ["nx","ny","cx","cy","pixel_size","detector_distance","wavelength","euler_angle_0","euler_angle_1","euler_angle_2"]
        if self._qmap_cache == {}:
            calculate = True
        else:
            for k in keys:
                exec "if self._qmap_cache[\"%s\"] != %s: calculate = True" % (k,k)
        if calculate:
            logger.info("Calculating qmap.")
            self._qmap_cache["qmap"] = generate_qmap(nx=nx, ny=ny, cx=cx, cy=cy, pixel_size=pixel_size, detector_distance=detector_distance, wavelength=wavelength,
                                                     euler_angle_0=euler_angle_0,euler_angle_1=euler_angle_1,euler_angle_2=euler_angle_2)
            for k in keys:
                exec "self._qmap_cache[\"%s\"] =  %s" % (k,k)
        return self._qmap_cache["qmap"]          

    #def get_map3d(self, map3d, dx, dx_req):
    #    map3d = None
    #    if dx > dx_req:
    #        logger.error("Finer real space sampling required for chosen geometry.")
    #        return
    #    # has map3d_fine the required real space grid?
    #    if map3d == None and abs(self.dX_fine/self.dX-1) < 0.001:
    #        # ok, we'll take the fine map
    #        map3d = self.map3d_fine
    #        logger.debug("Using the fine map for propagtion.")
    #        self._map3d = self.map3d_fine
    #    # do we have an interpolated map?
    #    if map3d == None and self._dX != None:
    #        # does it have the right spacing?
    #        if abs(self._dX/self.dX-1) < 0.001:
    #            # are the shapes of the original fine map and our current fine map the same?
    #            if numpy.all(numpy.array(self.map3d_fine.shape)==numpy.array(self._map3d_fine.shape)):
    #                # is the grid of the original fine map and the current fine map the same?
    #                if self.dX_fine == self._dX_fine:
    #                    # are the values of the original fine map and the cached fine map the same?
    #                    if numpy.all(self.map3d_fine==self._map3d_fine):
    #                        # ok, we take the cached map!
    #                        map3d = self._map3d
    #                        logger.debug("Using the cached interpolated map for propagtion.")
    #    # do we have to do interpolation?
    #    if map3d == None and self.dX_fine < self.dX:
    #        from scipy import ndimage
    #        f = self.dX_fine/self.dX
    #        N_mapfine = self.map3d_fine.shape[0]
    #        L_mapfine = (N_mapfine-1)*self.dX_fine
    #        N_map = int(numpy.floor((N_mapfine-1)*f))+1
    #        L_map = (N_map-1)*self.dX
    #        gt = numpy.float64(numpy.indices((N_map,N_map,N_map)))/float(N_map-1)*(N_mapfine-1)*L_map/L_mapfine
    #        map3d = ndimage.map_coordinates(self.map3d_fine, gt, order=3)
    #        # Cache interpolated data 
    #        self._map3d = map3d
    #        self._dX = N_mapfine/(1.*N_map)*self.dX_fine
    #        # Cace fine data for later decision whether or not the interpolated map can be used again
    #        self._map3d_fine = self.map3d_fine
    #        self._dX_fine = self.dX_fine
    #        logger.debug("Using a newly interpolated map for propagtion.")
    #    return map3d



def generate_absqmap(nx,ny,cx,cy,pixel_size,detector_distance,wavelength):
    X,Y = numpy.meshgrid(numpy.arange(nx),
                         numpy.arange(ny))
    # THIS CAST IS VERY IMPORTANT, in python A += B is not the same as A = A + B
    X = numpy.float64(X)
    Y = numpy.float64(Y)
    X -= cx
    Y -= cy
    return condortools.generate_absqmap(X,Y,pixel_size,detector_distance,wavelength)

def generate_qmap(nx,ny,cx,cy,pixel_size,detector_distance,wavelength,euler_angle_0=0.,euler_angle_1=0.,euler_angle_2=0.):
    X,Y = numpy.meshgrid(numpy.arange(nx),
                         numpy.arange(ny))
    # THIS CAST IS VERY IMPORTANT, in python A += B is not the same as A = A + B
    X = numpy.float64(X)
    Y = numpy.float64(Y)
    X -= cx
    Y -= cy
    return condortools.generate_qmap(X,Y,pixel_size,detector_distance,wavelength,euler_angle_0,euler_angle_1,euler_angle_2)
