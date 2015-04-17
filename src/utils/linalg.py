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

import sys, numpy, types, pickle, time, math
 
import logging
logger = logging.getLogger("Condor")
from log import log

def rotation(vector_or_matrix,E0,E1,E2):
    cos = numpy.cos
    sin = numpy.sin
    #Lsq = vector[0]**2+vector[1]**2+vector[2]**2
    M = numpy.array([[cos(E1)*cos(E2),
                      -cos(E0)*sin(E2)+sin(E0)*sin(E1)*cos(E2),
                      sin(E0)*sin(E2)+cos(E0)*sin(E1)*cos(E2)],
                     [cos(E1)*sin(E2),
                      cos(E0)*cos(E2)+sin(E0)*sin(E1)*sin(E2),
                      -sin(E0)*cos(E2)+cos(E0)*sin(E1)*sin(E2)],
                     [-sin(E1),
                      sin(E0)*cos(E1),
                      cos(E0)*cos(E1)]])
    rotated = numpy.dot(M,vector_or_matrix)
    #Lsq_rotated = rotated_vector[0]**2+rotated_vector[1]**2+rotated_vector[2]**2
    #if abs((Lsq - Lsq_rotated)/Lsq) > 0.001:
    #    print "ERROR: Rotation changes length!"
    #    print Lsq,Lsq_rotated
    return rotated

