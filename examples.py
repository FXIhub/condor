import os

this_dir = os.path.dirname(os.path.realpath(__file__))

def run_examples(on_travis=False):

    examples = [
        {
            "name": "PARTICLE IDEAL SPHERE (configfile)",
            "dir": this_dir + "/examples/configfile/particle_sphere",
            "cmd": "condor",
            "travis": True,
        },
        {
            "name": "PARTICLE IDEAL SPHEROID (configfile)",
            "dir": this_dir + "/examples/configfile/particle_spheroid",
            "cmd": "condor",
            "travis": True,
        },
        {
            "name": "PARTICLE MAP (configfile)",
            "dir": this_dir + "/examples/configfile/particle_map",
            "cmd": "condor",
            "travis": True,
        },
        {
            "name": "PARTICLE ATOMS (configfile)",
            "dir": this_dir + "/examples/configfile/particle_atoms",
            "cmd": "condor",
            "travis": True,
        },
        {
            "name": "PARTICLE CUSTOM MAP (script)",
            "dir": this_dir + "/examples/scripts/custom_map",
            "cmd": "python example.py",
            "travis": False,
        },
        {
            "name": "PARTICLE MAP EMD FETCH (script)",
            "dir": this_dir + "/examples/scripts/emd_fetch",
            "cmd": "python example.py",
            "travis": False,
        },
        {
            "name": "PARTICLE MAP MODELS (script)",
            "dir": this_dir + "/examples/scripts/particle_models",
            "cmd": "python example.py",
            "travis": True,
        },
        {
            "name": "PARTICLE ROTATIONS (script)",
            "dir": this_dir + "/examples/scripts/rotations",
            "cmd": "python example.py",
            "travis": True,
        },
        {
            "name": "PARTICLE SIMPLE (script)",
            "dir": this_dir + "/examples/scripts/simple",
            "cmd": "python example.py",
            "travis": True,
        },
        {
            "name": "PARTICLE ATOMS PDB FETCH (script)",
            "dir": this_dir + "/examples/scripts/pdb_fetch",
            "cmd": "python example.py",
            "travis": True,
        },
        {
            "name": "PARTICLE ATOMS PDB FILE (script)",
            "dir": this_dir + "/examples/scripts/pdb",
            "cmd": "python example.py",
            "travis": True,
        },
    ]
    
    if on_travis:
        examples = [e for e in examples if e["travis"]]

    nerrors = 0
        
    print "-"*100
    print ""
    for i,e in enumerate(examples):
        print ">>> Example %i/%i: %s" % (i+1, len(examples), e["name"])
        cmd = "cd %s; %s" % (e["dir"],e["cmd"])
        print cmd
        print ""
        error = os.system(cmd)
        print ""
        if error != 0:
            nerrors += 1
            raise Exception(">>> Example %i (%s) failed. Abort." % (i+1,e["name"]))
        else:
            print ">>> Success!"
        print ""
        print "-"*100
        print ""

    if nerrors == 0:
        print "SUCCESS: All examples finished successfully."
    else:
        print "ERROR: %i/%i example(s) failed." % (nerrors, len(examples))


if __name__ == "__main__":
    run_examples()
