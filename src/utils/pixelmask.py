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

class PixelMask:
    """
    CXI bitmasks follow the formalism described in Maia (2012) [#Maia2012]_ (see also `www.cxidb.org <http://www.cxidb.org>`_) and encode the status of detector pixels. In *Condor* these bitmasks are represented by arrays of unsigned 16 bit integer values. 

    The status of a pixel is a combination of pixel properties represented by the bits of the bitmask value. A pixel with no bit set (all bits are 0) represents an \"perfect\" pixel. No property is set.
    
    ``PIXEL_IS_PERFECT`` = 0
    
    Each bit that is set to 1 indicates that the respective property is assigned to that pixel.

    **Bit code:**

    ========================================  ============  =============================
    Name                                      Bit position  Bit in integer representation 
    ========================================  ============  =============================
    ``PIXEL_IS_INVALID``                      0th           0            
    ``PIXEL_IS_SATURATED``                    1st           2            
    ``PIXEL_IS_HOT``                          2nd           4            
    ``PIXEL_IS_DEAD``                         3rd           8 
    ``PIXEL_IS_SHADOWED``                     4th           16
    ``PIXEL_IS_IN_PEAKMASK``                  5th           32
    ``PIXEL_IS_TO_BE_IGNORED``                6th           64
    ``PIXEL_IS_BAD``                          7th           128
    ``PIXEL_IS_OUT_OF_RESOLUTION_LIMITS``     8th           256
    ``PIXEL_IS_MISSING``                      9th           512
    ``PIXEL_IS_NOISY``                        10th          1024
    ``PIXEL_IS_ARTIFACT_CORRECTED``           11th          2048
    ``PIXEL_FAILED_ARTIFACT_CORRECTION``      12th          4096
    ``PIXEL_IS_PEAK_FOR_HITFINDER``           13th          8192
    ``PIXEL_IS_PHOTON_BACKGROUND_CORRECTED``  14th          16384
    (not used)                                15th          32768
    ========================================  ============  =============================

    .. [#Maia2012] Filipe R. N. C. Maia, The Coherent X-ray Imaging Data Bank, Nature Methods 9, 854-855 (2012) `doi:10.1038/nmeth.2110 <http://dx.doi.org/doi:10.1038/nmeth.2110>`_.


    """
    # CXI pixelmask bits
    PIXEL_IS_PERFECT = 0
    PIXEL_IS_INVALID = 1
    PIXEL_IS_SATURATED = 2
    PIXEL_IS_HOT = 4
    PIXEL_IS_DEAD = 8
    PIXEL_IS_SHADOWED = 16
    PIXEL_IS_IN_PEAKMASK = 32
    PIXEL_IS_TO_BE_IGNORED = 64
    PIXEL_IS_BAD = 128
    PIXEL_IS_OUT_OF_RESOLUTION_LIMITS = 256
    PIXEL_IS_MISSING = 512
    PIXEL_IS_NOISY = 1024
    PIXEL_IS_ARTIFACT_CORRECTED = 2048
    PIXEL_FAILED_ARTIFACT_CORRECTION = 4096
    PIXEL_IS_PEAK_FOR_HITFINDER = 8192
    PIXEL_IS_PHOTON_BACKGROUND_CORRECTED = 16384
    # Cummulative bit options
    PIXEL_IS_IN_MASK = PIXEL_IS_INVALID |  PIXEL_IS_SATURATED | PIXEL_IS_HOT | PIXEL_IS_DEAD | PIXEL_IS_SHADOWED | PIXEL_IS_IN_PEAKMASK | PIXEL_IS_TO_BE_IGNORED | PIXEL_IS_BAD | PIXEL_IS_MISSING
