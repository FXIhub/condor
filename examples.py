import os

this_dir = os.path.dirname(os.path.realpath(__file__))

examples = [
    {
        "name": "PARTICLE IDEAL SPHERE (configfile)",
        "dir": this_dir + "/examples/configfile/particle_sphere",
        "cmd": "condor",
    },
    {
        "name": "PARTICLE IDEAL SPHEROID (configfile)",
        "dir": this_dir + "/examples/configfile/particle_spheroid",
        "cmd": "condor",
    },
    {
        "name": "PARTICLE MAP (configfile)",
        "dir": this_dir + "/examples/configfile/particle_map",
        "cmd": "condor",
    },
    {
        "name": "PARTICLE ATOMS (configfile)",
        "dir": this_dir + "/examples/configfile/particle_atoms",
        "cmd": "condor",
    },
#    {
#        "name": "PARTICLE CUSTOM MAP (script)",
#        "dir": this_dir + "/examples/scripts/custom_map",
#        "cmd": "python example.py",
#    },
#    {
#        "name": "PARTICLE MAP EMD FETCH (script)",
#        "dir": this_dir + "/examples/scripts/emd_fetch",
#        "cmd": "python example.py",
#    },
    {
        "name": "PARTICLE MAP MODELS (script)",
        "dir": this_dir + "/examples/scripts/particle_models",
        "cmd": "python example.py",
    },
    {
        "name": "PARTICLE ROTATIONS (script)",
        "dir": this_dir + "/examples/scripts/rotations",
        "cmd": "python example.py",
    },
    {
        "name": "PARTICLE SIMPLE (script)",
        "dir": this_dir + "/examples/scripts/simple",
        "cmd": "python example.py",
    },
    {
        "name": "PARTICLE ATOMS PDB FETCH (script)",
        "dir": this_dir + "/examples/scripts/pdb_fetch",
        "cmd": "python example.py",
   },
    {
        "name": "PARTICLE ATOMS PDB FILE (script)",
        "dir": this_dir + "/examples/scripts/pdb",
        "cmd": "python example.py",
   },

]

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
        raise Exception(">>> Example %i (%s) failed. Abort." % (i+1,e["name"]))
        nerrors += 1
    else:
        print ">>> Success!"
    print ""
    print "-"*100
    print ""

if nerrors == 0:
    print "SUCCESS: All examples finished successfully."
else:
    print "ERROR: %i/%i example(s) failed." % (nerrors, len(examples))
