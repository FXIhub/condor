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

class PixelMask:
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
