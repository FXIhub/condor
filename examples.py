import os

this_dir = os.path.dirname(os.path.realpath(__file__))

examples = [
    {
        "dir": this_dir + "/examples/configfile/particle_sphere",
        "cmd": "condor",
    },
    {
        "dir": this_dir + "/examples/configfile/particle_spheroid",
        "cmd": "condor",
    },
    {
        "dir": this_dir + "/examples/configfile/particle_map",
        "cmd": "condor",
    },
    {
        "dir": this_dir + "/examples/configfile/particle_molecule",
        "cmd": "condor",
    },
    {
        "dir": this_dir + "/examples/custom_map",
        "cmd": "condor",
    },
    {
        "dir": this_dir + "/examples/particle_models",
        "cmd": "python example.py",
    },
    {
        "dir": this_dir + "/examples/rotations",
        "cmd": "python example.py",
    },
    {
        "dir": this_dir + "/examples/simple",
        "cmd": "python example.py",
    },
]

    
for e in examples[4:]:
    cmd = "cd %s; %s" % (e["dir"],e["cmd"])
    print cmd
    os.system(cmd)
