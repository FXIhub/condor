import pylab, sys, numpy, types, pickle, time, math
import _config

PROPAGATION_MODE_PROJECTION = 0
PROPAGATION_MODE_3DSAMPLING = 1
PROPAGATOR_DIR = _config.PROPAGATOR_DIR

def init_configuration():
    # Load global dictionaries
    init_global_dictionaries()
    # Direct output to writable object
    commandline_out_deactivate()
    # Add path of propagator to sys.path
    sys.path.append(PROPAGATOR_DIR+"/utils")
    
def init_global_dictionaries():
    # Load scattering factors and atomic masses from Henke tables
    unpickle_scattering_factors()
    # Realative atomic compositions of certain material types (order: H,C,N.O,P,S)
    global DICT_atomic_composition
    DICT_atomic_composition = {'protein':[86,52,13,15,0,3],
                               'cell':[23,3,1,10,0,1], # Bergh et al. 2008
                               'latex':[1,1,0,0,0,0], 
                               'water':[2,0,0,1,0,0], 
                               'dna':[11,10,4,6,1,0],
                               'lipid':[69,36,0,6,1,0],
                               'genophore':[205,134,38,48,3,6],
                               'virus':[72.43,49.85,16.32,24.49,2.57,1.39],
                               'mimivirus':[23,3,1,10,0,1],
                               'carboxysome':[0.51,0.30,0.07,0.10,0.,0.02],
                               'sucrose':[22,12,0,11,0,0]}
    # Estimated mass densities of certain material types
    global DICT_massdensity
    DICT_massdensity = {'protein':1350,
                        'cell':1000,
                        'latex':1050,
                        'water':998,
                        'Au':19300,
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

    # Physical constants [SI-units]
    global DICT_physical_constants
    DICT_physical_constants = {'e':1.602176487E-19,
                               'c':299792458,
                               'h':6.62606896E-34,
                               're':2.8179402894E-15,
                               'barn':1E-28,
                               'u':1.66053886E-27}


def unpickle_scattering_factors():
    global DICT_atomic_mass
    DICT_atomic_mass = {}
    global DICT_scattering_factors
    DICT_scattering_factors = {}
    ELEMENTS_FILE = open('%s/elements.dat' % PROPAGATOR_DIR,'r')
    DICT_atomic_mass,DICT_scattering_factors = pickle.load(ELEMENTS_FILE)
    F_MIN_ENERGY_EV = 0
    F_MAX_ENERGY_EV = 0
    for var in DICT_scattering_factors.values():
        if F_MIN_ENERGY_EV < var[0,0] or F_MIN_ENERGY_EV == 0: F_MIN_ENERGY_EV = var[0,0]
        if F_MAX_ENERGY_EV > var[-1,0] or F_MAX_ENERGY_EV == 0: F_MAX_ENERGY_EV = var[-1,0]

class _WritableObject:
    """ Class that can be assigned to stdout in order to switch off output to commandline and instead to write it to self.content"""
    def __init__(self):
        self.content = []
    def write(self, string):
        self.content.append(string)

def commandline_out_deactivate():
    global OUT
    OUT = _WritableObject()

def commandline_out_activate():
    global OUT
    OUT = sys.stdout

