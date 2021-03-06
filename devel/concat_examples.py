#!/usr/bin/env python
"""
Create example configuration files by concatenation
"""
from __future__ import print_function, absolute_import # Compatibility with python 2 and 3
import os,sys

repodir = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..')

def concatenate_files(infilenames,outfilename,extra_newlines=0):
    lines = []
    for infilename in infilenames:
        with open(infilename, "r") as f:
            lines += f.readlines()
        if extra_newlines > 0:
            lines += ["\n"]*extra_newlines
    with open(outfilename, "w") as f:
        f.writelines(lines)

if __name__ == "__main__":

    print("Concatenating example configuration files...")

    src = os.path.join(repodir, 'examples', 'configfile', 'source.conf')
    det = os.path.join(repodir, 'examples', 'configfile', 'detector.conf')

    for m in ["particle_sphere","particle_spheroid","particle_map","particle_atoms"]:
        par = os.path.join(repodir, 'examples', 'configfile', '%s.conf' % m)
        infilenames = [src, par, det]
        d = os.path.join(repodir, 'examples', 'configfile', '%s' % m)
        if not os.path.exists(d):
            os.mkdir(d)
        outfilename = os.path.join(d, "condor.conf")
        print("Concatenating " + outfilename + "...")
        concatenate_files(infilenames, outfilename, 2)
    
    print("Done.")
