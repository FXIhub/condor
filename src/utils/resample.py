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

import sys, numpy

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
