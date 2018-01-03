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
"""
Reading, writing and converting configuration in different representations
"""

from __future__ import print_function, absolute_import # Compatibility with python 2 and 3
import os, numpy, tempfile, copy
try:
    import ConfigParser as configparser
except ImportError:
    # In Python 3, configparser needs to be installed with pip
    import configparser

import logging
logger = logging.getLogger(__name__)

from .log import log_and_raise_error,log_warning,log_info,log_debug
import condor

def read_configfile(configfile):
    """
    Read configuration file to dictionary
    """
    config = configparser.ConfigParser()
    with open(configfile,"r") as f:
        config.readfp(f)
        confDict = {}
        for section in config.sections(): 
            confDict[section] = {}
            c = config.items(section)
            for (key,value) in c:
                confDict[section][key] = _estimate_class(value)
    return confDict

def write_configfile(configdict, filename):
    """
    Write configuration file from a dictionary
    """
    ls = ["# Configuration file\n# Automatically written by Configuration instance\n\n"]
    for section_name,section in configdict.items():
        if isinstance(section,dict):
            ls.append("[%s]\n" % section_name)
            for variable_name,variable in section.items():
                if (hasattr(variable, '__len__') and (not isinstance(variable, str))) or isinstance(variable, list):
                    ls.append("%s=%s\n" % (variable_name,_list_to_str(variable)))
                else:
                    ls.append("%s=%s\n" % (variable_name,str(variable)))
            ls.append("\n")
    with open(filename, "w") as f:
        f.writelines(ls)

def read_configdict(configdict):
    C = {}
    for k,v in configdict.items():
        if isinstance(v, dict):
            v_new = read_configdict(v)
        else:
            v_new = _estimate_class(v)
        C[k] = v_new
    return C

def _estimate_class(var):
    v = _estimate_type(var)
    if isinstance(v,str):
        v = v.replace(" ","")
        if v.startswith("[") and v.endswith("]"):
            v = _str_to_list(v)
            for i in range(len(v)):
                v[i] = os.path.expandvars(v[i]) if isinstance(v[i], str) else v[i] 
        elif v.startswith("{") and v.endswith("}"):
            v = v[1:-1].split(",")
            v = [w for w in v if len(w) > 0]
            d = {}
            for w in v:
                key,value = w.split(":")
                value = _estimate_type(value)
                d[key] = value
            v = d
        else:
            if v.startswith("$"):
                v = os.path.expandvars(v)
    return v
        
def _estimate_type(var):
    if not isinstance(var, str):
        return var
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

def _str_to_list(s):
    if s.startswith("[") and s.endswith("]"):
        if s[1:-1].startswith("[") and s[1:-1].endswith("]"):
            return _str_to_list(s[1:-1])
        else:
            l = s[1:-1].split(",")
            l = [_estimate_type(w) for w in l if len(w) > 0]
            return l
    else:
        return s
       
def _list_to_str(L):
    if (hasattr(L, '__len__') and (not isinstance(L, str))) or isinstance(L, list):
        s = ""
        for l in L:
            s += _list_to_str(l)
            s += ","
        s = "[" + s[:-1] + "]"
        return s
    else:
        return str(L)
            
def _conf_to_spsim_opts(D_source,D_particle,D_detector,ndim=2,qn=None,qmax=None):
    if ndim == 2:
        if qn is not None or qmax is not None:
            log_warning(logger, "As ndim=2 the passed values for qn and qmax take no effect.")
    if ndim == 3:
        if qn is None and qmax is None:
            log_and_raise_error(logger, "As ndim=3 both qn and qmax must be not None.")
            return
    import spsim
    # Create temporary file for pdb file
    tmpf_pdb = tempfile.NamedTemporaryFile(mode='w+', suffix='.conf', prefix='tmp_spsim', dir=None, delete=False)
    tmpf_pdb_name = tmpf_pdb.name
    tmpf_pdb.close()
    # Write pdb file
    mol = spsim.get_molecule_from_atoms(D_particle["atomic_numbers"], D_particle["atomic_positions"])
    spsim.write_pdb_from_mol(tmpf_pdb_name, mol)
    spsim.free_mol(mol)
    # Start with default spsim configuration
    opts = spsim.set_defaults()
    # Create temporary file for spsim configuration
    tmpf_conf = tempfile.NamedTemporaryFile(mode='w+', suffix='.conf', prefix='tmp_spsim', dir=None, delete=False)
    # Write string sequence from configuration dicts
    s = []
    s += "# THIS FILE WAS CREATED AUTOMATICALLY BY CONDOR\n"
    s += "# Temporary configuration file for spsim\n"
    s += "verbosity_level = 0;\n"
    s += "number_of_dimensions = %i;\n" % ndim
    s += "number_of_patterns = 1;\n"
    s += "origin_to_com = 1;\n"
    s += "input_type = \"pdb\";\n"
    #s += "pdb_filename = \"%s\";\n" % D_particle["pdb_filename"]
    s += "pdb_filename = \"%s\";\n" % tmpf_pdb_name
    if ndim == 2:
        D = D_detector["distance"]
        Lx = D_detector["pixel_size"] * D_detector["nx"]
        Ly = D_detector["pixel_size"] * D_detector["ny"]
    else:
        k0 = 2. * numpy.pi / D_source["wavelength"]
        D = qn / 2. * D_detector["pixel_size"] * k0 / qmax
        Lx = Ly = Lz = D_detector["pixel_size"] * qn
    s += "detector_distance = %.12e;\n" % D
    s += "detector_width = %.12e;\n" % Lx 
    s += "detector_height = %.12e;\n" % Ly
    if ndim == 3:
        s += "detector_depth = %.12e;\n" % Lz        
    s += "detector_pixel_width = %.12e;\n" % D_detector["pixel_size"]
    s += "detector_pixel_height = %.12e;\n" % D_detector["pixel_size"]
    if ndim == 3:
        s += "detector_pixel_depth = %.12e;\n" % D_detector["pixel_size"]
    if ndim == 2:
        s += "detector_center_x = %.12e;\n" % (D_detector["pixel_size"] * (D_detector["cx"] - (D_detector["nx"]-1)/2.))
        s += "detector_center_y = %.12e;\n" % (D_detector["pixel_size"] * (D_detector["cy"] - (D_detector["ny"]-1)/2.))
    else:
        s += "detector_center_x = 0;\n"
        s += "detector_center_y = 0;\n"
        s += "detector_center_z = 0;\n"
    s += "detector_binning = 1;\n"
    s += "experiment_wavelength = %.12e;\n" % D_source["wavelength"]
    s += "experiment_beam_intensity = %.12e;\n" % D_particle["intensity"]
    s += "experiment_polarization = \"ignore\";\n" # polarization correction will be done in Condor if needed (see experiment.py)
    #s += "use_cuda = 0;\n"
    intrinsic_rotation = condor.utils.rotation.Rotation(values=D_particle["extrinsic_quaternion"],formalism="quaternion")
    intrinsic_rotation.invert()
    e0, e1, e2 = intrinsic_rotation.get_as_euler_angles("zxz")
    if not numpy.isfinite(e0):
        print("ERROR: phi is not finite")
    if not numpy.isfinite(e1):
        print("ERROR: theta is not finite")
    if not numpy.isfinite(e2):
        print("ERROR: psi is not finite")
    s += "phi = %.12e;\n" % e0
    s += "theta = %.12e;\n" % e1
    s += "psi = %.12e;\n" % e2
    s += "random_orientation = 0;\n"
    # Write string sequence to file
    tmpf_conf.writelines(s)
    # Close temporary file
    tmpf_conf_name = tmpf_conf.name
    tmpf_conf.close()
    # Read configuration into options struct
    spsim.read_options_file(tmpf_conf_name, opts)
    # This deletes the temporary files
    os.unlink(tmpf_pdb_name)
    os.unlink(tmpf_conf_name)
    return opts

