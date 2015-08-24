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

import os, glob, numpy, pickle, sys, re

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
        
    with open('%s/sf.dat' % (outpath),'w') as f:
        pickle.dump(S,f)

def pickle_atomic_standard_weights_and_numbers(inpath, outpath):
    with open(inpath + "/standard_weights.txt", "r") as f:
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

    with open(outpath + "/sw.dat", "w") as f:
        pickle.dump(W,f)

    with open(outpath + "/z.dat", "w") as f:
        pickle.dump(Z,f)
