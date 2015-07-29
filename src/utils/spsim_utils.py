import os, tempfile

import spsim

def conf_to_opts(D_source,D_particle,D_detector):
    # Start with default spsim configuration
    opts = spsim.set_defaults()
    # Create temporary file for spsim configuration
    tmpf = tempfile.NamedTemporaryFile(mode='w+b', bufsize=-1, suffix='.conf', prefix='tmp_spsim', dir=None, delete=False)
    # Write string sequence from configuration dicts
    s = []
    s += "# THIS FILE WAS CREATED AUTOMATICALLY BY CONDOR\n"
    s += "# Temporary configuration file for spsim\n"
    s += "verbosity_level = 0;\n"
    s += "number_of_dimensions = 2;\n"
    s += "number_of_patterns = 1;\n"
    s += "input_type = \"pdb\";\n"
    s += "pdb_filename = \"%s\";\n" % D_particle["pdb_filename"]
    s += "detector_distance = %.6e;\n" % D_detector["distance"]
    s += "detector_width = %.6e;\n" % (D_detector["pixel_size"] * D_detector["nx"]) 
    s += "detector_height = %.6e;\n" % (D_detector["pixel_size"] * D_detector["ny"])
    s += "detector_pixel_width = %.6e;\n" % D_detector["pixel_size"]
    s += "detector_pixel_height = %.6e;\n" % D_detector["pixel_size"]
    s += "detector_center_x = %.6e;\n" % (D_detector["pixel_size"] * (D_detector["cx"] - (D_detector["nx"]-1)/2.))
    s += "detector_center_y = %.6e;\n" % (D_detector["pixel_size"] * (D_detector["cy"] - (D_detector["ny"]-1)/2.))
    s += "detector_binning = 1;\n"
    s += "experiment_wavelength = %.6e;\n" % D_source["wavelength"]
    s += "experiment_beam_intensity = %.6e;\n" % D_particle["intensity"]
    #s += "phi = %.6e;\n" % D_particle["euler_angle_0"]
    #s += "theta = %.6e;\n" % D_particle["euler_angle_1"]
    #s += "psi = %.6e;\n" % D_particle["euler_angle_2"]
    s += "phi = 1.0;\n"
    s += "theta = 1.0;\n"
    s += "psi = 1.0;\n"
    s += "random_orientation = 0;\n"
    # Write string sequence to file
    tmpf.writelines(spsim_conf)
    # Close temporary file
    tmpf_name = tmpf.name
    tmpf.close()
    # Read configuration into options struct
    spsim.read_options_file(tmpf_name, opts)
    # This deletes the temporary file
    os.unlink(tmpf_name)
    return opts

def opts_to_mol(opts):
    mol = spsim.get_molecule(opts)

def atoms_to_mol(atomic_number = None, atomic_position = None):
    mol = spsim.alloc_mol()
    for j,(pos0,pos1,pos2) in zip(atomic_number, atomic_position.reshape(len(atomic_number), 3)):
        spsim.add_atom_to_mol(mol, int(j), pos0, pos1, pos2)
    return mol

def mol_to_atoms(mol):
    pos_img = spsim.sp_image_alloc(mol.natoms, 3, 1)
    spsim.array_to_image(mol.pos, pos_img)
    pos = pos_img.image.real[:,:].copy()
    spsim.sp_image_free(pos_img)
    anum_img = spsim.sp_image_alloc(mol.natoms, 1, 1)
    spsim.iarray_to_image(mol.atomic_number, anum_img)
    anum = numpy.int32(anum_img.image.real[:,:].copy())
    spsim.sp_image_free(anum_img)
    return [anum,pos]
