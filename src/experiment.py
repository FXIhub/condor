# -----------------------------------------------------------------------------------------------------
# CONDOR
# Simulator for diffractive single-particle imaging experiments with X-ray lasers
# http://xfel.icm.uu.se/condor/
# -----------------------------------------------------------------------------------------------------
# Copyright 2014 Max Hantke, Filipe R.N.C. Maia, Tomas Ekeberg
# Condor is distributed under the terms of the GNU General Public License
# -----------------------------------------------------------------------------------------------------
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but without any warranty; without even the implied warranty of
# merchantability or fitness for a pariticular purpose. See the
# GNU General Public License for more details.
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
# -----------------------------------------------------------------------------------------------------
# General note:
# All variables are in SI units by default. Exceptions explicit by variable name.
# -----------------------------------------------------------------------------------------------------


# TO DO:
# Take into account illumination profile

import numpy, os, sys, copy

import logging
logger = logging.getLogger(__name__)

import sample
from condor.utils.log import log,log_execution_time
from condor.utils.log import log_and_raise_error,log_warning,log_info,log_debug
import condor.utils.config
from condor.utils.pixelmask import PixelMask
import condor.utils.sphere_diffraction
import condor.utils.spheroid_diffraction
import condor.utils.scattering_vector
import condor.utils.resample
from condor.utils.rotation import Rotation
import condor.sample
import condor.particle

class Experiment:
    def __init__(self, source, sample, detector):
        self.source   = source
        self.sample   = sample
        self.detector = detector
        self._qmap_cache = {}

    def get_conf(self):
        conf = {}
        conf.update(self.source.get_conf())
        conf.update(self.sample.get_conf())
        conf.update(self.detector.get_conf())
        return conf

    @log_execution_time(logger)
    def propagate(self, save_map3d = False, save_qmap = False):
        N = self.source.number_of_shots
        O = {"source":{}, "sample":{}, "detector":{}}
        O_particles = {}
        N_particles_max = 0

        for i in range(N):
            log_debug(logger, "Calculation diffraction pattern (%i/%i). (PROGSTAT)" % (i+1, N))

            Os = self._propagate_single(save_map3d=save_map3d, save_qmap=save_qmap)

            N_particles = len(Os["sample"]["particles"])
            log_info(logger, "%i/%i (%i particle%s)" % (i+1, N, N_particles, "s" if N_particles > 1 else ""))

            for k in [k for k in Os.keys() if k not in ["source", "sample", "detector"]]:
                if k not in O:
                    O[k] = [None for _ in range(N)]
                O[k][i] = Os[k]
            # Source
            for k,v in Os["source"].items():
                if k not in O["source"]:
                    O["source"][k] = [None for _ in range(N)]
                O["source"][k][i] = Os["source"][k]
            # Sample
            for k in [k for k in Os["sample"].keys() if k != "particles"]:
                if k not in O["sample"]:
                    O["sample"][k] = [None for _ in range(N)]
                O["sample"][k][i] = Os["sample"][k]
            for j,p in zip(range(len(Os["sample"]["particles"])), Os["sample"]["particles"]):
                for k,v in p.items():
                    if k not in O_particles:
                        O_particles[k] = [[] for _ in range(N)]
                    O_particles[k][i].append(v)
            N_particles_max = max([N_particles_max, N_particles])
            # Detector
            for k,v in Os["detector"].items():
                if k not in O["detector"]:
                    O["detector"][k] = [None] * N
                O["detector"][k][i] = Os["detector"][k]
        # Make arrays for particle parameters
        for k,v in O_particles.items():
            s = [N, N_particles_max]
            v0 = numpy.array(v[0])
            s_item = list(v0.shape)
            s_item.pop(0)
            s = s + s_item
            A = numpy.zeros(shape=s, dtype=v0.dtype)
            for i in range(len(v)):
                vi = numpy.array(v[i])
                d = len(vi.shape)
                if d == 1:
                    A[i,:vi.shape[0]] = vi[:]
                elif d == 2:
                    A[i,:vi.shape[0],:] = vi[:,:]
                elif d == 3:
                    A[i,:vi.shape[0],:,:] = vi[:,:,:]
                elif d == 4:
                    A[i,:vi.shape[0],:,:,:] = vi[:,:,:,:]
            O["sample"][k] = A
        for k in O["source"].keys(): O["source"][k] = numpy.array(O["source"][k])
        for k in O["detector"].keys(): O["detector"][k] = numpy.array(O["detector"][k])
        for k in [k for k in O.keys() if k not in ["source", "sample", "detector"]]:
            O[k] = numpy.array(O[k])
        return O

    @log_execution_time(logger)
    def _propagate_single(self, save_map3d = False, save_qmap = False):
        
        # Iterate objects
        D_source   = self.source.get_next()
        D_sample   = self.sample.get_next()
        D_detector = self.detector.get_next()

        # Pull out variables
        nx                  = D_detector["nx"]
        ny                  = D_detector["ny"]
        cx                  = D_detector["cx"]
        cy                  = D_detector["cy"]
        pixel_size          = D_detector["pixel_size"]
        detector_distance   = D_detector["distance"]
        Omega_p             = D_detector["solid_angle_pixel"]
        wavelength          = D_source["wavelength"]

        F_singles    = []
        qmap_singles = []
        # Calculate patterns of all single particles individually
        for D_particle in D_sample["particles"]:
            p  = D_particle["_class_instance"]
            # Intensity at interaction point
            pos  = D_particle["position"]
            D_particle["intensity"] = self.source.get_intensity(pos, "ph/m2", pulse_energy=D_source["pulse_energy"])
            I_0 = D_particle["intensity"]
            # Calculate primary wave amplitude
            # F0 = sqrt(I_0 Omega_p) 2pi/wavelength^2
            F0 = numpy.sqrt(I_0*Omega_p)*2*numpy.pi/wavelength**2
            D_particle["F0"] = F0
            # 3D Orientation
            extrinsic_rotation = Rotation(values=D_particle["extrinsic_quaternion"], formalism="quaternion")
            
            # UNIFORM SPHERE
            if isinstance(p, condor.particle.ParticleSphere):
                # Refractive index
                dn = p.material.get_dn(wavelength)
                # Scattering vectors
                qmap = self.get_qmap(nx=nx, ny=ny, cx=cx, cy=cy, pixel_size=pixel_size, detector_distance=detector_distance, wavelength=wavelength, extrinsic_rotation=None)
                q = numpy.sqrt(qmap[:,:,1]**2+qmap[:,:,2]**2)
                # Intensity scaling factor
                R = D_particle["diameter"]/2.
                V = 4/3.*numpy.pi*R**3
                K = (F0*V*abs(dn))**2
                # Pattern
                F = condor.utils.sphere_diffraction.F_sphere_diffraction(K, q, R)

            # UNIFORM SPHEROID
            elif isinstance(p, condor.particle.ParticleSpheroid):
                # Refractive index
                dn = p.material.get_dn(wavelength)
                # Scattering vectors
                qmap = self.get_qmap(nx=nx, ny=ny, cx=cx, cy=cy, pixel_size=pixel_size, detector_distance=detector_distance, wavelength=wavelength, extrinsic_rotation=None, order="xyz")
                qx = qmap[:,:,0]
                qy = qmap[:,:,1]
                # Intensity scaling factor
                R = D_particle["diameter"]/2.
                V = 4/3.*numpy.pi*R**3
                K = (F0*V*abs(dn))**2
                a = condor.utils.spheroid_diffraction.to_spheroid_semi_diameter_a(D_particle["diameter"], D_particle["flattening"])
                c = condor.utils.spheroid_diffraction.to_spheroid_semi_diameter_c(D_particle["diameter"], D_particle["flattening"])
                # Pattern
                # Spheroid axis before rotation
                v0 = numpy.array([0.,1.,0.])
                v1 = extrinsic_rotation.rotate_vector(v0)
                theta = numpy.arcsin(v1[2])
                phi   = numpy.arctan2(-v1[0],v1[1])
                F = condor.utils.spheroid_diffraction.F_spheroid_diffraction(K, qx, qy, a, c, theta, phi)

            # MAP
            elif isinstance(p, condor.particle.ParticleMap):
                # Refractive index
                dn = p.material.get_dn(wavelength)
                # Scattering vectors (the nfft requires order z,y,x)
                qmap = self.get_qmap(nx=nx, ny=ny, cx=cx, cy=cy, pixel_size=pixel_size, detector_distance=detector_distance, wavelength=wavelength, extrinsic_rotation=extrinsic_rotation, order="zyx")
                # Intensity scaling factor
                R = D_particle["diameter"]/2.
                V = 4/3.*numpy.pi*R**3
                # Generate map
                # We multiply by 0.99 to prevent numerical issues
                dx_required  = self.detector.get_resolution_element_r(wavelength, cx=cx, cy=cy, center_variation=False) * 0.99
                dx_suggested = self.detector.get_resolution_element_r(wavelength, center_variation=True) * 0.99
                map3d, dx = p.get_map3d(D_particle, dx_required, dx_suggested)
                log_debug(logger, "Sampling of map: dx_required = %e m, dx_suggested = %e m, dx = %e m" % (dx_required, dx_suggested, dx))
                map3d_dn = numpy.array(map3d, dtype=numpy.complex128) * dn
                if save_map3d:
                    D_particle["map3d_dn"] = map3d_dn
                    D_particle["dx"] = dx
                # Rescale and shape qmap for nfft
                qmap_scaled = dx * qmap / (2. * numpy.pi)
                qmap_shaped = qmap_scaled.reshape(qmap_scaled.shape[0]*qmap_scaled.shape[1], 3)
                # For Jing:
                #D_particle["qmap3d"] = detector.generate_qmap_ori(nfft_scaled=True)
                # Check inputs
                invalid_mask = ((qmap_shaped>=-0.5) * (qmap_shaped<0.5)) == False
                if (invalid_mask).sum() > 0:
                    qmap_shaped[invalid_mask] = 0.
                    log_warning(logger, "%i invalid pixel positions." % invalid_mask.sum())
                log_debug(logger, "Map3d input shape: (%i,%i,%i), number of dimensions: %i, sum %f" % (map3d_dn.shape[0], map3d_dn.shape[1], map3d_dn.shape[2], len(list(map3d_dn.shape)), abs(map3d_dn).sum()))
                if (numpy.isfinite(abs(map3d_dn))==False).sum() > 0:
                    log_warning(logger, "There are infinite values in the dn map of the object.")
                log_debug(logger, "Scattering vectors shape: (%i,%i); Number of dimensions: %i" % (qmap_shaped.shape[0], qmap_shaped.shape[1], len(list(qmap_shaped.shape))))
                if (numpy.isfinite(qmap_shaped)==False).sum() > 0:
                    log_warning(logger, "There are infinite values in the scattering vectors.")
                # NFFT
                fourier_pattern = log_execution_time(logger)(condor.utils.nfft.nfft)(map3d_dn, qmap_shaped)
                # Check output - masking in case of invalid values
                if (invalid_mask).sum() > 0:
                    fourier_pattern[numpy.any(invalid_mask)] = numpy.nan
                # reshaping
                fourier_pattern = numpy.reshape(fourier_pattern, (qmap_scaled.shape[0], qmap_scaled.shape[1]))
                log_debug(logger, "Got pattern of %i x %i pixels." % (fourier_pattern.shape[1], fourier_pattern.shape[0]))
                F = F0 * fourier_pattern * dx**3

            # MOLECULE
            elif isinstance(p, condor.particle.ParticleMolecule):
                # Import only here (otherwise errors if spsim library not installed)
                import spsim
                # Create options struct
                opts = condor.utils.config.conf_to_opts(D_source, D_particle, D_detector)
                spsim.write_options_file("./spsim.confout",opts)
                # Create molecule struct
                mol = spsim.get_molecule_from_atoms(D_particle["atomic_numbers"], D_particle["atomic_positions"])
                spsim.write_pdb_from_mol("./mol.pdbout", mol)
                # Calculate diffraction pattern
                pat = spsim.simulate_shot(mol, opts)
                # Extract complex Fourier values from spsim output
                F_img = spsim.make_cimage(pat.F, pat.rot, opts)
                phot_img = spsim.make_image(opts.detector.photons_per_pixel, pat.rot, opts)
                F = numpy.sqrt(abs(phot_img.image[:])) * numpy.exp(1.j * numpy.angle(F_img.image[:]))
                spsim.sp_image_free(F_img)
                spsim.sp_image_free(phot_img)
                # Extract qmap from spsim output
                qmap_img = spsim.sp_image_alloc(3,D_detector["nx"], D_detector["ny"])
                spsim.array_to_image(pat.HKL_list, qmap_img)
                qmap = numpy.zeros(shape=(D_detector["ny"], D_detector["nx"], 3))
                qmap[:,:,0] = qmap_img.image.real[:,:,0]
                qmap[:,:,1] = qmap_img.image.real[:,:,1]
                qmap[:,:,2] = qmap_img.image.real[:,:,2]
                spsim.sp_image_free(qmap_img)
                spsim.free_diffraction_pattern(pat)
                spsim.free_output_in_options(opts)
                
            else:
                log_and_raise_error(logger, "No valid particles initialized.")
                sys.exit(0)

            if save_qmap:
                qmap_singles.append(qmap)

            F_singles.append(F)

        F_tot = numpy.zeros_like(F)
        # Superimpose patterns
        for D_particle,F in zip(D_sample["particles"], F_singles):
            v = D_particle["position"]
            F_tot = F_tot + F * numpy.exp(-1.j*(v[0]*qmap[:,:,0]+v[1]*qmap[:,:,1]+v[2]*qmap[:,:,2])) 
        I_tot, M_tot, IXxX_tot, MXxX_tot = self.detector.detect_photons(abs(F_tot)**2)
        if self.detector.binning is not None:
            FXxX_tot, MXxX_tot = condor.utils.resample.downsample(F_tot, self.detector.binning, mode="integrate", 
                                                                  mask2d0=M_tot, bad_bits=PixelMask.PIXEL_IS_IN_MASK, min_N_pixels=1)
        M_tot_binary = M_tot == 0
        MXxX_tot_binary = None if MXxX_tot is None else (MXxX_tot == 0)
        
        O = {}
        O["source"]            = D_source
        O["sample"]            = D_sample
        O["detector"]          = D_detector

        O["fourier_pattern"]   = F_tot
        O["intensity_pattern"] = I_tot
        O["mask_binary"]       = M_tot_binary
        O["mask"]              = M_tot

        #O["qmap"] = qmap_singles
        
        if self.detector.binning is not None:
            O["fourier_pattern_xxx"]   = FXxX_tot
            O["intensity_pattern_xxx"] = IXxX_tot
            O["mask_xxx"]              = MXxX_tot
            O["mask_xxx_binary"]       = MXxX_tot_binary
            
        return O

    @log_execution_time(logger)
    def get_qmap(self, nx, ny, cx, cy, pixel_size, detector_distance, wavelength, extrinsic_rotation=None, order="xyz"):
        calculate = False
        if self._qmap_cache == {}:
            calculate = True
        else:
            calculate = calculate or nx != self._qmap_cache["nx"]
            calculate = calculate or ny != self._qmap_cache["ny"]
            calculate = calculate or cx != self._qmap_cache["cx"]
            calculate = calculate or cy != self._qmap_cache["cy"]
            calculate = calculate or pixel_size != self._qmap_cache["pixel_size"]
            calculate = calculate or detector_distance != self._qmap_cache["detector_distance"]
            calculate = calculate or wavelength != self._qmap_cache["wavelength"]
            calculate = calculate or not extrinsic_rotation.similar(self._qmap_cache["extrinsic_rotation"])
        if calculate:
            log_debug(logger,  "Calculating qmap")
            Y,X = numpy.meshgrid(numpy.arange(ny), numpy.arange(nx), indexing="ij")
            X = numpy.float64(X)
            Y = numpy.float64(Y)
            X -= cx
            Y -= cy
            self._qmap_cache = {
                "qmap"              : condor.utils.scattering_vector.generate_qmap(X, Y, pixel_size, detector_distance, wavelength, extrinsic_rotation=extrinsic_rotation, order=order),
                "nx"                : nx,
                "ny"                : ny,
                "cx"                : cx,
                "cy"                : cy,
                "pixel_size"        : pixel_size,
                "detector_distance" : detector_distance,
                "wavelength"        : wavelength,
                "extrinsic_rotation": copy.deepcopy(extrinsic_rotation),
            }            
        return self._qmap_cache["qmap"]          

    def get_resolution(self, wavelength = None, cx = None, cy = None, pos="corner", convention="full_period"):
        if wavelength is None:
            wavelength = self.source.photon.get_wavelength()
        dx = self.detector.get_resolution_element(wavelength, cx=cx, cy=cy, pos=pos)
        if convention == "full_period":
            return dx*2
        elif convention == "half_period":
            return dx
        else:
            log_and_raise_error(logger, "Invalid input: convention=%s. Must be either \"full_period\" or \"half_period\"." % convention)
            return

    def get_linear_sampling_ratio(self, wavelength = None, particle_diameter = None):
        """
        Returns the linear sampling ratio :math:`o` of the diffraction pattern:

        | :math:`o=\\frac{D\\lambda}{dp}` 

        | :math:`D`: Detector distance
        | :math:`p`: Detector pixel size (edge length)
        | :math:`\\lambda`: Photon wavelength 
        | :math:`d`: Particle diameter

        """
        if wavelength is None:
            wavelength = self.source.photon.get_wavelength()
        detector_distance = self.detector.distance
        if particle_diameter is None:
            pm = self.sample.get_particle_models()
            N = len(pm.keys()) 
            if N == 0:
                log_and_raise_error(logger, "You need to specify a particle_diameter because no particle model is defined.")
                return
            elif N > 1:
                log_and_raise_error(logger, "The particle_diameter is ambiguous because more than one particle model is defined.")
                return
            else:
                particle_diameter = pm.values()[0]
        pN = utils.diffraction.nyquist_pixel_size(wavelength, detector_distance, particle_diameter)
        pD = self.detector.pixel_size
        ratio = pN/pD
        return ratio
        
    def get_fresnel_number(self, wavelength):
        pass
                           
    
    # ------------------------------------------------------------------------------------------------
    # Caching of map3d might be interesting to implement again in the future
    #def get_map3d(self, map3d, dx, dx_req):
    #    map3d = None
    #    if dx > dx_req:
    #        condor.CONDOR_logger.error("Finer real space sampling required for chosen geometry.")
    #        return
    #    # has map3d_fine the required real space grid?
    #    if map3d == None and abs(self.dX_fine/self.dX-1) < 0.001:
    #        # ok, we'll take the fine map
    #        map3d = self.map3d_fine
    #        condor.CONDOR_logger.debug("Using the fine map for proagtion.")
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
    #                        condor.CONDOR_logger.debug("Using the cached interpolated map for propagtion.")
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
    #        condor.CONDOR_logger.debug("Using a newly interpolated map for propagtion.")
    #    return map3d
    # ------------------------------------------------------------------------------------------------
