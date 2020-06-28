#!/usr/bin/env python
from __future__ import print_function, absolute_import # Compatibility with python 2 and 3
import os

repodir = os.path.dirname(os.path.realpath(__file__))
examplesdir_configfile = os.path.join(repodir, "examples", "configfile")
examplesdir_scripts = os.path.join(repodir, "examples", "scripts")
examplesdir_publication = os.path.join(repodir, "examples_publication")

def run_examples(on_travis=False):

    examples = [
        {
            "name": "PARTICLE IDEAL SPHERE (configfile)",
            "dir": os.path.join(examplesdir_configfile, "particle_sphere"),
            "cmd": "rm -f condor.cxi; condor",
            "travis": True,
        },
        {
            "name": "PARTICLE IDEAL SPHEROID (configfile)",
            "dir": os.path.join(examplesdir_configfile, "particle_spheroid"),
            "cmd": "rm -f condor.cxi; condor",
            "travis": True,
        },
        {
            "name": "PARTICLE MAP (configfile)",
            "dir": os.path.join(examplesdir_configfile, "particle_map"),
            "cmd": "rm -f condor.cxi; condor",
            "travis": True,
        },
        {
            "name": "PARTICLE ATOMS (configfile)",
            "dir": os.path.join(examplesdir_configfile, "particle_atoms"),
            "cmd": "rm -f condor.cxi; condor",
            "travis": True,
        },
        {
            "name": "PARTICLE CUSTOM MAP (script)",
            "dir": os.path.join(examplesdir_scripts, "custom_map"),
            "cmd": "rm -f condor.cxi; python example.py",
            "travis": False,
        },
        {
            "name": "PARTICLE MAP EMD FETCH (script)",
            "dir": os.path.join(examplesdir_scripts, "emd_fetch"),
            "cmd": "rm -f condor.cxi; python example.py",
            "travis": False,
        },
        {
            "name": "PARTICLE MAP MODELS (script)",
            "dir": os.path.join(examplesdir_scripts, "particle_models"),
            "cmd": "rm -f condor.cxi; python example.py",
            "travis": True,
        },
        {
            "name": "PARTICLE ROTATIONS (script)",
            "dir": os.path.join(examplesdir_scripts, "rotations"),
            "cmd": "rm -f condor.cxi; python example.py",
            "travis": True,
        },
        {
            "name": "PARTICLE SIMPLE (script)",
            "dir": os.path.join(examplesdir_scripts, "simple"),
            "cmd": "rm -f condor.cxi; python example.py",
            "travis": True,
        },
        {
            "name": "PARTICLE ATOMS PDB FETCH (script)",
            "dir": os.path.join(examplesdir_scripts, "pdb_fetch"),
            "cmd": "rm -f condor.cxi; python example.py",
            "travis": True,
        },
        {
            "name": "PARTICLE ATOMS PDB FILE (script)",
            "dir": os.path.join(examplesdir_scripts, "pdb"),
            "cmd": "rm -f condor.cxi; python example.py",
            "travis": True,
        },
        {
            "name": "PARTICLE MAP MULTIPLE MATERIALS (script)",
            "dir": os.path.join(examplesdir_scripts, "multiple_materials"),
            "cmd": "python example.py",
            "travis": False,
        },
        {
            "name": "PARTICLE MAP DIFFRACTION SPACE 3D INTERPOLATION (script)",
            "dir": os.path.join(examplesdir_scripts, "diffraction_3d"),
            "cmd": "python example.py",
            "travis": False,
        },
        {
            "name": "DIFFRACTION SPACE 3D SIMULATION (script)",
            "dir": os.path.join(examplesdir_scripts, "full_fourier_volume"),
            "cmd": "python example.py",
            "travis": True,
        },
        {
            "name": "PUBLICATION EXAMPLE A: PARTICLE ATOMS GROEL (configfile)",
            "dir": os.path.join(examplesdir_publication, "a"),
            "cmd": "rm -f condor.cxi; condor",
            "travis": True,
        },
        {
            "name": "PUBLICATION EXAMPLE B: PARTICLE MAP EMD1144 (configfile)",
            "dir": os.path.join(examplesdir_publication, "b"),
            "cmd": "rm -f condor.cxi; condor",
            "travis": False,
        },
        {
            "name": "PUBLICATION EXAMPLE B: PARTICLE MAP EMD1144 (script)",
            "dir": os.path.join(examplesdir_scripts, "publication_example_b"),
            "cmd": "rm -f condor.cxi; python example.py",
            "travis": False,
        },        
    ]
    
    if on_travis:
        examples = [e for e in examples if e["travis"]]

    nerrors = 0
        
    print("-"*100)
    print("")
    for i,e in enumerate(examples):
        print(">>> Example %i/%i: %s" % (i+1, len(examples), e["name"]))
        cmd = "cd %s; %s" % (e["dir"],e["cmd"])
        print(cmd)
        print("[start output]")
        error = os.system(cmd)
        print("[end output]")
        if error != 0:
            nerrors += 1
            raise Exception(">>> Example %i (%s) failed. Abort." % (i+1,e["name"]))
        else:
            print(">>> Success!")
        print("")
        print("-"*100)
        print("")

    if nerrors == 0:
        print("SUCCESS: All examples finished successfully.")
    else:
        print("ERROR: %i/%i example(s) failed." % (nerrors, len(examples)))


if __name__ == "__main__":
    run_examples()
