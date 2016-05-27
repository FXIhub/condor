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

import urllib2
import StringIO
import gzip

import numpy

import scipy.ndimage
import scipy.stats
import scipy.ndimage.measurements


import logging
logger = logging.getLogger(__name__)

import log

def fetch_map(emd_id):
    url = "ftp://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-%s/map/emd_%s.map.gz" % (str(emd_id),str(emd_id))
    log.log_debug(logger, "Downloading file for EMDID %s from URL %s" % (emd_id, url))
    filename = "./emd_%s.map" % str(emd_id)
    response = urllib2.urlopen(url)
    compressedFile = StringIO.StringIO()
    compressedFile.write(response.read())
    compressedFile.seek(0)
    decompressedFile = gzip.GzipFile(fileobj=compressedFile, mode='rb')
    log.log_debug(logger, "Download of %s ended." % (filename))
    with open(filename, 'w') as outfile:
        outfile.write(decompressedFile.read())
    return read_map(filename)
    
def read_map(filename):
    log.log_warning(logger, "Automatic scaling of EM maps may not be reliable. Please make sure to check your map after using this functionality.")
    # CCP4 map file format
    # http://www.ccp4.ac.uk/html/maplib.html
    with open(filename, "rb") as f:
        # 1024 bytes header
        header_buf = f.read(1024)
        temp_int32 = numpy.frombuffer(header_buf, dtype="int32")
        temp_float32 = numpy.frombuffer(header_buf, dtype="float32")
        #1      NC              # of Columns    (fastest changing in map)
        #2      NR              # of Rows
        #3      NS              # of Sections   (slowest changing in map)
        NC = temp_int32[0]
        NR = temp_int32[1]
        NS = temp_int32[2]
        if NC != NR or NR != NS:
            log.log_and_raise_error(logger, "Cannot read a map with unequal dimensions")
        N = NC
        #4      MODE            Data type
        #                  0 = envelope stored as signed bytes (from
        #                      -128 lowest to 127 highest)
        #                  1 = Image     stored as Integer*2
        #                  2 = Image     stored as Reals
        #                  3 = Transform stored as Complex Integer*2
        #                  4 = Transform stored as Complex Reals
        #                  5 == 0	
        #
        #                  Note: Mode 2 is the normal mode used in
        #                        the CCP4 programs. Other modes than 2 and 0
        #                        may NOT WORK        
        MODE = temp_int32[3]
        dtype = ["int8", "int16", "float32", None, "complex64", "int8"][MODE]
        if MODE == 3:
            log.log_and_raise_error(logger, "Map file data type \"MODE=%i\" is not implemented yet." % MODE)
        if MODE not in [0,2,5]:
            log.log_warning(logger, "Map file data type \"MODE=%i\" not supported yet and may not work reliably." % MODE)
        #11      X length        Cell Dimensions (Angstroms)
        #12      Y length                     "
        #13      Z length                     "
        dX = temp_float32[10]/float(N)*1E-10
        dY = temp_float32[11]/float(N)*1E-10
        dZ = temp_float32[12]/float(N)*1E-10
        if dX != dY or dY != dZ:
            log.log_and_raise_error(logger, "Cannot read a map with unequal voxel dimensions")
        #17      MAPC            Which axis corresponds to Cols.  (1,2,3 for X,Y,Z)
        #18      MAPR            Which axis corresponds to Rows   (1,2,3 for X,Y,Z)
        #19      MAPS            Which axis corresponds to Sects. (1,2,3 for X,Y,Z)
        MAPC = temp_int32[16]
        MAPR = temp_int32[17]
        MAPS = temp_int32[18]
        #24      NSYMBT          Number of bytes used for storing symmetry operators
        NSYMBT = temp_int32[23]
        if NSYMBT > 0:
            log.log_and_raise_error(logger, "Omitting symmetry operations in map file.")
            f.read(NSYMBT)
        # The remaining bytes are data
        raw_data = f.read()
        raw_data = numpy.frombuffer(raw_data, dtype=dtype)
        # Now we need to project onto the right Z-Y-X array grid
        S,R,C = numpy.meshgrid(numpy.arange(NS), numpy.arange(NR), numpy.arange(NC), indexing='ij')
        S = S.flatten()
        R = R.flatten()
        C = C.flatten()
        if MAPC == 1:
            X = C
            Xlen = NC
        elif MAPC == 2:
            Y = C
            Ylen = NC
        elif MAPC == 3:
            Z = C
            Zlen = NC
        if MAPR == 1:
            X = R
            Xlen = NR
        elif MAPR == 2:
            Y = R
            Ylen = NR
        elif MAPR == 3:
            Z = R
            Zlen = NR
        if MAPS == 1:
            X = S
            Xlen = NS
        elif MAPS == 2:
            Y = S
            Ylen = NS
        elif MAPS == 3:
            Z = S
            Zlen = NS
        i = Z*(Ylen*Xlen) + Y*(Xlen) + X
        i.sort()
        data = numpy.zeros(Zlen*Ylen*Xlen, dtype=dtype)
        data[:] = raw_data[i]
        data = data.reshape((Zlen,Ylen,Xlen))
    return data, dX


def preproc_map_auto(map3d_raw, ed_water, ed_particle, water_layer=0.1):
    N = map3d_raw.shape[0]
    Nr = (N-1)/2.
    Nw = int(numpy.round(Nr*water_layer))
    Nw = 1 if Nw==0 else Nw
    Z,Y,X = numpy.meshgrid(numpy.arange(N)-Nr,
                           numpy.arange(N)-Nr,
                           numpy.arange(N)-Nr, indexing="ij")
    R = numpy.sqrt(X**2+Y**2+Z**2)
    mask = R < Nr
    water_shell = (abs(R-(R*mask).max()) <= Nw)*mask
    S0 = numpy.sort(map3d_raw[mask])
    c = S0.min() + (S0.max()-S0.min())/2.
    S1 = S0[S0<c]
    S2 = S0[S0>=c]
    v_water = scipy.stats.mode(S1)[0][0]
    v_particle_min = scipy.stats.mode(S2)[0][0]
    threshold = v_water + (v_particle_min - v_water)*0.1
    map3d_bin = mask * map3d_raw > threshold
    labeled_array, num_features = scipy.ndimage.measurements.label(map3d_bin)
    map3d_bin2 = labeled_array == labeled_array[N/2,N/2,N/2]
    map3d_bin_filled = scipy.ndimage.morphology.binary_fill_holes(map3d_bin2)
    v_particle = numpy.mean(map3d_raw[map3d_bin_filled])
    ed_map3d = (ed_water + (map3d_raw-v_water)/(v_particle-v_water)*(ed_particle-ed_water))/ed_particle * map3d_bin_filled

    #import h5py
    #with h5py.File("/Users/hantke/test.h5","w") as f:
    #    f["mask"] = mask
    #    f["ed_map3d"] = ed_map3d
    #    f["map3d_raw"] = map3d_raw
    #    f["map3d_bin_filled"] = map3d_bin_filled
    return ed_map3d

def preproc_map_manual(map3d_raw, offset, factor):
    map3d = (map3d_raw + offset) * factor
