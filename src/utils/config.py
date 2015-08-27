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

import os, numpy, ConfigParser, tempfile, copy

import logging
logger = logging.getLogger(__name__)

from log import log_and_raise_error,log_warning,log_info,log_debug
import condor

def load_config(config = None, default = None):
    if config is None:
        config = {}
    if default is None:
        default = {}
    if isinstance(config,str):
        confDict = read_configfile(config)
    else:
        confDict = config
        
    if isinstance(default,str):
        defDict = read_configfile(default)
    else:
        defDict = default
            
    # set unspecified to default
    for sectionName in defDict.keys():
        if sectionName not in confDict.keys():
            self.confDict[sectionName] = {}
            log_info(logger, "Add section %s with as section does not exist." % (sectionName))
        for variableName in [n for n in defDict[sectionName].keys() if n not in confDict[sectionName].keys()]:
            confDict[sectionName][variableName] = defDict[sectionName][variableName]
            log_info(logger, "Add variable %s with default value %s to configuration section %s as variable does not exist." % (variableName,str(defDict[sectionName][variableName]),sectionName))

    # set condor variables
    for sectionName,section in confDict.items():
        for variableName,variable in section.items():
            if isinstance(variable, basestring):
                if variable.startswith("CONDOR_"):
                    exec("confDict[sectionName][variableName] = condor.%s" % variable)

    return confDict

def write_config(confDict, filename):
    ls = ["# Configuration file\n# Automatically written by Configuration instance\n\n"]
    for section_name,section in confDict.items():
        if isinstance(section,dict):
            ls.append("[%s]\n" % section_name)
            for variable_name,variable in section.items():
                if (hasattr(variable, '__len__') and (not isinstance(variable, str))) or isinstance(variable, list):
                    ls.append("%s=%s\n" % (variable_name,list_to_str(variable)))
                else:
                    ls.append("%s=%s\n" % (variable_name,str(variable)))
            ls.append("\n")
    with open(filename, "w") as f:
        f.writelines(ls)
 
def estimate_type(var):
    #first test bools
    if var.lower() == 'true':
        return True
    elif var.lower() == 'false':
        return False
    elif var.lower() == 'none':
        return None
    else:
        #int
        try:
            return int(var)
        except ValueError:
            pass
        #float
        try:
            return float(var)
        except ValueError:
            pass
        #string
        try:
            return str(var)
        except ValueError:
            raise NameError('Something messed up autocasting var %s (%s)' % (var, type(var)))

def read_configfile(configfile):
    config = ConfigParser.ConfigParser()
    with open(configfile,"r") as f:
        config.readfp(f)
        confDict = {}
        for section in config.sections(): 
            confDict[section] = {}
            c = config.items(section)
            for (key,value) in c:
                v = estimate_type(value)
                if isinstance(v,str):
                    v = v.replace(" ","")
                    if v.startswith("[") and v.endswith("]"):
                        v = str_to_list(v)
                        for i in range(len(v)):
                            v[i] = os.path.expandvars(v[i]) if isinstance(v[i], str) else v[i] 
                    elif v.startswith("{") and v.endswith("}"):
                        v = v[1:-1].split(",")
                        v = [w for w in v if len(w) > 0]
                        d = {}
                        for w in v:
                            key,value = w.split(":")
                            value = estimate_type(value)
                            if value.startswith("$"):
                                value = os.path.expandvars(value)
                            d[key] = value
                        v = d
                    else:
                        if v.startswith("$"):
                            v = os.path.expandvars(v)
                confDict[section][key] = v
    return confDict

def str_to_list(s):
    if s.startswith("[") and s.endswith("]"):
        if s[1:-1].startswith("[") and s[1:-1].endswith("]"):
            return str_to_list(s[1:-1])
        else:
            l = s[1:-1].split(",")
            l = [estimate_type(w) for w in l if len(w) > 0]
            return l
    else:
        return s
       
def list_to_str(L):
    if (hasattr(L, '__len__') and (not isinstance(L, str))) or isinstance(L, list):
        s = ""
        for l in L:
            s += list_to_str(l)
            s += ","
        s = "[" + s[:-1] + "]"
        return s
    else:
        return str(L)
            
# SPSIM
def conf_to_opts(D_source,D_particle,D_detector):
    import spsim
    # Create temporary file for pdb file
    tmpf_pdb = tempfile.NamedTemporaryFile(mode='w+b', bufsize=-1, suffix='.conf', prefix='tmp_spsim', dir=None, delete=False)
    tmpf_pdb_name = tmpf_pdb.name
    tmpf_pdb.close()
    # Write pdb file
    mol = spsim.get_molecule_from_atoms(D_particle["atomic_numbers"], D_particle["atomic_positions"])
    spsim.write_pdb_from_mol(tmpf_pdb_name, mol)
    spsim.free_mol(mol)
    # Start with default spsim configuration
    opts = spsim.set_defaults()
    # Create temporary file for spsim configuration
    tmpf_conf = tempfile.NamedTemporaryFile(mode='w+b', bufsize=-1, suffix='.conf', prefix='tmp_spsim', dir=None, delete=False)
    # Write string sequence from configuration dicts
    s = []
    s += "# THIS FILE WAS CREATED AUTOMATICALLY BY CONDOR\n"
    s += "# Temporary configuration file for spsim\n"
    s += "verbosity_level = 0;\n"
    s += "number_of_dimensions = 2;\n"
    s += "number_of_patterns = 1;\n"
    s += "input_type = \"pdb\";\n"
    #s += "pdb_filename = \"%s\";\n" % D_particle["pdb_filename"]
    s += "pdb_filename = \"%s\";\n" % tmpf_pdb_name
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
    intrinsic_rotation = condor.utils.rotation.Rotation(values=D_particle["extrinsic_quaternion"],formalism="quaternion")
    intrinsic_rotation.invert()
    e0, e1, e2 = intrinsic_rotation.get_as_euler_angles("zxz")
    s += "phi = %.6e;\n" % e0
    s += "theta = %.6e;\n" % e1
    s += "psi = %.6e;\n" % e2
    s += "random_orientation = 0;\n"
    # Write string sequence to file
    tmpf_conf.writelines(s)
    # Close temporary file
    tmpf_conf_name = tmpf_conf.name
    tmpf_conf.close()
    # Read configuration into options struct
    spsim.read_options_file(tmpf_conf_name, opts)
    # This deletes the temporary files

    os.system("cp %s /Users/hantke/foo.pdb" % tmpf_pdb_name)
    os.unlink(tmpf_pdb_name)
    os.unlink(tmpf_conf_name)
    return opts


