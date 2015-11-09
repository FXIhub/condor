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

import urllib2
import StringIO
import gzip

import numpy

import logging
logger = logging.getLogger(__name__)

import log

def fetch_map(emd_id):
    url = "ftp://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-%s/map/emd_%s.map.gz" % (str(emd_id),str(emd_id))
    filename = "./emd_%s.map" % str(emd_id)
    response = urllib2.urlopen(url)
    compressedFile = StringIO.StringIO()
    compressedFile.write(response.read())
    compressedFile.seek(0)
    decompressedFile = gzip.GzipFile(fileobj=compressedFile, mode='rb')
    with open(filename, 'w') as outfile:
        outfile.write(decompressedFile.read())
    return read_map(filename)
    
def read_map(filename):
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
        dtype = ["int8", "uint16", "float32", None, "complex64", "int8"][MODE]
        if MODE == 3:
            log.log_and_raise_error(logger, "Map file data type \"MODE=%i\" is not implemented yet." % MODE)
        if MODE not in [0,2,5]:
            log.log_warning(logger, "Map file data type \"MODE=%i\" not supported yet and may not work." % MODE)
        #11      X length        Cell Dimensions (Angstroms)
        #12      Y length                     "
        #13      Z length                     "
        dX = temp_float32[10]*1E-9
        dY = temp_float32[11]*1E-9
        dZ = temp_float32[12]*1E-9
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


