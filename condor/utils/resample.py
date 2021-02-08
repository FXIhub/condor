# -----------------------------------------------------------------------------------------------------
# CONDOR
# Simulator for diffractive single-particle imaging experiments with X-ray lasers
# http://xfel.icm.uu.se/condor/
# -----------------------------------------------------------------------------------------------------
# Copyright 2016 Max Hantke, Filipe R.N.C. Maia, Tomas Ekeberg
# Condor is distributed under the terms of the BSD 2-Clause License
# -----------------------------------------------------------------------------------------------------
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# 
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# -----------------------------------------------------------------------------------------------------
# General note:
# All variables are in SI units by default. Exceptions explicit by variable name.
# -----------------------------------------------------------------------------------------------------

from __future__ import print_function, absolute_import # Compatibility with python 2 and 3
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
        print("ERROR: %s is not a valid mode." % mode)
        return
    factor = int(round(factor0))
    if factor == 1:
        if mask2d0 is None:
            return array2d0
        else:
            return [array2d0,mask2d0]
    array2d = numpy.array(array2d0,dtype=array2d0.dtype)
    if mask2d0 is None:
        mask2d = None
    else:
        mask2d = numpy.array(mask2d0,dtype="int16")
    Nx = array2d.shape[1]
    Ny = array2d.shape[0]
    if mode == "pick": 
        Y,X = numpy.indices(array2d0.shape)
        pick = ((Y%factor == 0)*(X%factor == 0))
        print(pick.shape)
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
        A = A.flatten()
        Y,X = numpy.indices((Ny,Nx))
        Y = Y.flatten()
        X = X.flatten()
        Y //= factor
        X //= factor
        superp = Y*Nx_new+X
        superp_order = superp.argsort()
        A = A[superp_order]
        A = A.reshape((Nx_new*Ny_new,factor*factor))
        if mask2d is None:
            B = A.sum(1)
            return B.reshape((Ny_new,Nx_new))
        if mask2d is not None:
            AM = numpy.zeros(shape=(Ny,Nx),dtype="int16")
            AM[:mask2d.shape[0],:mask2d.shape[1]] = mask2d[:,:]
            AM = AM.flatten()
            AM = AM[superp_order]
            AM = AM.reshape((Nx_new*Ny_new,factor*factor))
            if bad_bits is None:
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
            B[BN >= min_N_pixels] = B[BN >= min_N_pixels] * factor*factor /numpy.float64(BN[BN >= min_N_pixels])
            return [B.reshape((Ny_new,Nx_new)),BM.reshape((Ny_new,Nx_new))]
