#!/usr/bin/env python
"""
Pickle constants data form text tables
"""
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
import os, glob, numpy, pickle, sys, re
import codecs

here = os.path.dirname(os.path.realpath(__file__))


def pickle_atomic_scattering_factors(inpath, outpath): 
    S = {}
    for infile in glob.glob(os.path.join("./", '%s/*.nff' % (inpath)) ):
        regexp = re.search("%s/([a-z]+).nff$" % (inpath),infile)
        el = regexp.group(1).capitalize()
        S[el] = []
        with open(infile,"r") as f:
            lines = f.readlines()
            lines = lines[1:]
            for line in lines:
                arg = line.split()
                if float(arg[1]) > 0:
                    S[el].append([float(arg[0]),float(arg[1]),float(arg[2])])
        S[el] = numpy.array(S[el])
        
    with open('%s/sf.dat' % (outpath),'wb') as f:
        pickle.dump(S,f)

def pickle_atomic_standard_weights_and_numbers(inpath, outpath):
    with codecs.open(inpath + "/standard_weights.txt", "r") as f:
        lines = f.readlines()
    # Skip the header
    lines = [l for l in lines if not l.startswith("#")]
    # Read atom numbers and standard weights
    W = {}
    Z = {}
    for l in lines:
        strs = l.split("\t")
        z = int(strs[0])
        s = strs[1]
        w = strs[3]
        w = w.replace("(","")
        w = w.replace(")","")
        w = w.replace("[","")
        w = w.replace("]","")
        w = w.replace("#",".")
        w = float(w)
        W[s] = w
        Z[s] = z

    with open(outpath + "/sw.dat", "wb") as f:
        pickle.dump(W,f)

    with open(outpath + "/z.dat", "wb") as f:
        pickle.dump(Z,f)



if __name__ == "__main__":
    
    # Atomic scattering factors from the Henke tables
    # B.L. Henke, E.M. Gullikson, and J.C. Davis. X-ray interactions: photoabsorption, scattering, transmission, and reflection at E=50-30000 eV, Z=1-92
    # Atomic Data and Nuclear Data Tables Vol. 54 (no.2), 181-342 (July 1993).
    # http://henke.lbl.gov/optical_constants/asf.html
    print('Generate data file with atomic scattering constants...')
    pickle_atomic_scattering_factors(here + "/sf", here)
    print('Done.')
    
    # Standard atomic weights from the IUPAC tables
    # Atomic weights of the elements 2011 (IUPAC Technical Report) Michael E. Wieser et al., Pure and Applied Chemistry. Volume 85, Issue 5, Pages 1047-1078
    # ISSN (Online) 1365-3075, ISSN (Print) 0033-4545
    # DOI: 10.1351/PAC-REP-13-03-02, April 2013
    # Data loaded from: http://www.chem.qmul.ac.uk/iupac/AtWt/ table 2 (2015/07/01)
    print('Generate data file with atomic standard weight constants...')
    pickle_atomic_standard_weights_and_numbers(here + "/sw", here)
    print('Done.')

