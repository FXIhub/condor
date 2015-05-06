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

import sys, os, numpy, types, pickle, time, math, logging, ConfigParser
logging.basicConfig(format='%(levelname)s: %(message)s')
logger = logging.getLogger('Condor')

def init_configuration():
    # Some global configuration variables
    #=============
    # Load global dictionaries
    init_global_dictionaries()
    # Add path of Condor to sys.path
    #global CONDOR_DIR
    CONDOR_DIR = os.path.dirname(os.path.realpath(__file__))
    sys.path.append(CONDOR_DIR+"/utils")
    global _CONDOR_DEFAULT_PDB
    _CONDOR_DEFAULT_PDB = CONDOR_DIR + "/data/DNA.pdb"
    
def init_global_dictionaries():
    # Load scattering factors and atomic masses from Henke tables
    unpickle_scattering_factors()
    # Realative atomic compositions of certain material types (order: H,C,N,O,P,S,Au)
    global DICT_atomic_composition
    DICT_atomic_composition = {'protein':[86,52,13,15,0,3,0],
                               'cell':[23,3,1,10,0,1,0], # Bergh et al. 2008
                               'latex':[1,1,0,0,0,0,0], 
                               'water':[2,0,0,1,0,0,0], 
                               'dna':[11,10,4,6,1,0,0],
                               'lipid':[69,36,0,6,1,0,0],
                               'genophore':[205,134,38,48,3,6,0],
                               'virus':[72.43,49.85,16.32,24.49,2.57,1.39,0],
                               'carboxysome':[0.51,0.30,0.07,0.10,0.,0.02,0],
                               'sucrose':[22,12,0,11,0,0,0],
                               'gold':[0,0,0,0,0,0,1]}
    # Estimated mass densities of certain material types
    global DICT_massdensity
    DICT_massdensity = {'protein':1350,
                        'cell':1000,
                        'latex':1050,
                        'water':998,
                        'gold':19300,
                        'dna':1700,
                        'lipid':1000,
                        'genophore':1560,
                        'virus':1381,
                        'mimivirus':1100,
                        'carboxysome':1250,
                        'sucrose':1587}
    # More documentation needed!
    # The following material types should be defined more properly:
    # - 'virus': density = 1455 (Filipe's webpage)
    # - 'carboxysome' density = 1250 (guessed by Dirk, atomic composition deduced assuming protein and water being the only two components)
    global DICT_atomic_number
    DICT_atomic_number = {'H':1,
                          'He':2,
                          'Li':3,
                          'Be':4,
                          'B':5,
                          'C':6,
                          'N':7,
                          'O':8,
                          'F':9,
                          'Ne':10,
                          'Na':11,
                          'Mg':12,
                          'Al':13,
                          'Si':14,
                          'P':15,
                          'S':16,
                          'Cl':17,
                          'Ar':18,
                          'K':19,
                          'Ca':20,
                          'Sc':21,
                          'Ti':22,
                          'V':23,
                          'Cr':24,
                          'Mn':25,
                          'Fe':26,
                          'Co':27,
                          'Ni':28,
                          'Cu':29,
                          'Zn':30,
                          'Ga':31,
                          'Ge':32,
                          'As':33,
                          'Se':34,
                          'Br':35,
                          'Kr':36,
                          'Rb':37,
                          'Sr':38,
                          'Y':39,
                          'Zr':40,
                          'Nb':41,
                          'Mo':42,
                          'Tc':43,
                          'Ru':44,
                          'Rh':45,
                          'Pd':46,
                          'Ag':47,
                          'Cd':48,
                          'In':49,
                          'Sn':50,
                          'Sb':51,
                          'Te':52,
                          'I':53,
                          'Xe':54,
                          'Cs':55,
                          'Ba':56,
                          'La':57,
                          'Ce':58,
                          'Pr':59,
                          'Nd':60,
                          'Pm':61,
                          'Sm':62,
                          'Eu':63,
                          'Gd':64,
                          'Tb':65,
                          'Dy':66,
                          'Ho':67,
                          'Er':68,
                          'Tm':69,
                          'Yb':70,
                          'Lu':71,
                          'Hf':72,
                          'Ta':73,
                          'W':74,
                          'Re':75,
                          'Os':76,
                          'Ir':77,
                          'Pt':78,
                          'Au':79,
                          'Hg':80,
                          'Tl':81,
                          'Pb':82,
                          'Bi':83,
                          'Po':84,
                          'At':85,
                          'Rn':86,
                          'Fr':87,
                          'Ra':88,
                          'Ac':89,
                          'Th':90,
                          'Pa':91,
                          'U':92,
                          'Np':93,
                          'Pu':94,
                          'Am':95,
                          'Cm':96,
                          'Bk':97,
                          'Cf':98,
                          'Es':99,
                          'Fm':100,
                          'Md':101,
                          'No':102,
                          'Lr':103,
                          'Rf':104,
                          'Db':105,
                          'Sg':106,
                          'Bh':107,
                          'Hs':108,
                          'Mt':109,
                          'Ds':110,
                          'Rg':111,
                          'Cp':112,
                          'Uut':113,
                          'Uuq':114,
                          'Uup':115,
                          'Uuh':116,
                          'Uus':117,
                          'Uuo':118}

def unpickle_scattering_factors():
    global DICT_atomic_mass
    DICT_atomic_mass = {}
    global DICT_scattering_factors
    DICT_scattering_factors = {}
    this_dir = os.path.dirname(os.path.realpath(__file__))
    ELEMENTS_FILE = open('%s/data/elements.dat' % this_dir,'r')
    DICT_atomic_mass,DICT_scattering_factors = pickle.load(ELEMENTS_FILE)
    F_MIN_ENERGY_EV = 0
    F_MAX_ENERGY_EV = 0
    for var in DICT_scattering_factors.values():
        if F_MIN_ENERGY_EV < var[0,0] or F_MIN_ENERGY_EV == 0: F_MIN_ENERGY_EV = var[0,0]
        if F_MAX_ENERGY_EV > var[-1,0] or F_MAX_ENERGY_EV == 0: F_MAX_ENERGY_EV = var[-1,0]


def check_input(keys,req_keys,opt_keys,verbose=False):
    missing_keys = [k for k in req_keys if not isinstance(k,list) and (k not in keys)]
    # Required keys from different alternative combination sets (list in req_keys)
    for l in [l for k in keys if isinstance(k,list)]:
        tmp_miss = []
        # Alternatives (non-exclusive)
        for a in l:
            if isinstance(a,list):
                # Combined requirement
                com = [c for c in a if c not in keys]
                if len(com) > 0:
                    tmp_miss.append(com)
            else:
                # Single requirement
                if a not in keys:
                    tmp_miss.append(a)
        if len(tmp_miss) == len(l):
            missing_keys.append(tmp_miss)
    # Required keys from the optional keys (keys that are needed in combined sets, indicated by tuples in opt_keys list)
    req_opt_keys = []
    for t in [t for t in opt_keys if isinstance(t,tuple)]:
        if list(set(t).intersection(keys)):
            for tt in t:
                if tt not in req_opt_keys:
                    req_opt_keys.append(tt)
    all_keys = req_keys + opt_keys
    def sublist_elements(l):
        l_new = []
        for k0 in l:
            if isinstance(k0,list):
                for k1 in k0:
                    l_new.append(k1)
            else:
                l_new.append(k0)
        return l_new
    illegal_keys = [k for k in keys if k not in all_keys]
    illegal_keys = [k for k in illegal_keys if k not in sublist_elements(all_keys)]
    illegal_keys = [k for k in illegal_keys if k not in sublist_elements(sublist_elements(all_keys))]
    if verbose:
        for illegal_key in illegal_keys:
            print "Illegal key: %s" % illegal_key
        if len(missing_keys) > 0:
            print "Missing key(s):"
        for missing_key in missing_keys:
            if isinstance(missing_key,list):
                print " - Alternatives:"
                for missing_key_alternative in missing_key:
                    if isinstance(missing_key_alternative,list):
                        s = "  - ["
                        for missing_key_alternative_component in missing_key_alternative:
                            s += missing_key_alternative_component + " + "
                        s += "]"
                        print s
                    else:
                        print ("  - " + missing_key_alternative)
            else:
                print (" - " + missing_key)
    return missing_keys,illegal_keys

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
                    if "[" == v[0] and "]" == v[-1]:
                        v = v[1:-1].split(",")
                        v = [w for w in v if len(w) > 0]
                        for i in range(len(v)):
                            if '$' in v[i]:
                                v[i] = os.path.expandvars(v[i])
                            v[i] = estimate_type(v[i]) 
                    else:
                        if '$' in v:
                            v = os.path.expandvars(v)
                confDict[section][key] = v
    return confDict

class Configuration:
    def __init__(self,config={},default={},verbose=False):
        self.verbose = verbose
        if isinstance(config,str):
            self.confDict = read_configfile(config)
        else:
            self.confDict = config
        
        if isinstance(default,str):
            defDict = read_configfile(default)
        else:
            defDict = default
            
        self.set_unspecified_to_default(defDict)
        self.replace_condor_variables_by_values()
        
    def set_unspecified_to_default(self,defaultDict):
        for sectionName in defaultDict.keys():
            if sectionName not in self.confDict.keys():
                self.confDict[sectionName] = {}
                logger.info("Add section %s with as section does not exist." % (sectionName))
            for variableName in [n for n in defaultDict[sectionName].keys() if n not in self.confDict[sectionName].keys()]:
                self.confDict[sectionName][variableName] = defaultDict[sectionName][variableName]
                logger.info("Add variable %s with default value %s to configuration section %s as variable does not exist." % (variableName,str(defaultDict[sectionName][variableName]),sectionName))

    def replace_condor_variables_by_values(self):
        for sectionName,section in self.confDict.items():
            for variableName,variable in section.items():
                if isinstance(variable, basestring):
                    if variable[:len("_CONDOR_")] == "_CONDOR_":
                        exec("self.confDict[sectionName][variableName] = %s" % variable)

    def write_to_file(self,filename):
        ls = ["# Configuration file\n# Automatically written by Configuration instance\n\n"]
        for section_name,section in self.confDict.items():
            if isinstance(section,dict):
                ls.append("[%s]\n" % section_name)
                for variable_name,variable in section.items():
                    if (hasattr(variable, '__len__') and (not isinstance(variable, str))) or isinstance(variable, list):
                        ls.append("%s=%s\n" % (variable_name,list_to_str(variable)))
                    else:
                        ls.append("%s=%s\n" % (variable_name,str(variable)))
                ls.append("\n")
        s = open(filename,"w")
        s.writelines(ls)
        s.close()        

        
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
            
