# Import library packages
#------------
import pylab, sys, numpy, types, pickle, time, math
SPOWPY_PATH = "/home/hantke/programs/propagator/"
IMGTOOLS_PATH = "/home/hantke/pythonscripts/tools"
sys.path.append(IMGTOOLS_PATH)

def unpickle_scattering_factors():

    global DICT_atomic_mass
    global DICT_scattering_factors
    ELEMENTS_FILE = open('%s/elements.dat' % SPOWPY_PATH,'r')
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


# Define global dictionaries
#---------------------------

# Typical realative atomic compositions (order: H,C,N.O,P,S), collected from Filipe Maia's webpage and paper Bergh et al. 2008
DICT_atomic_composition = {'protein':[86,52,13,15,0,3],'cell':[23,3,1,10,0,1],'latex':[1,1,0,0,0,0],'water':[2,0,0,1,0,0],'dna':[11,10,4,6,1,0],'lipid':[69,36,0,6,1,0],'genophore':[205,134,38,48,3,6],'virus':[72.43,49.85,16.32,24.49,2.57,1.39]}

# Typical realative atomic compositions (order: H,C,N.O,P,S)
DICT_massdensity = {'protein':1350,'cell':1000,'latex':1050,'water':998,'Au':19300,'dna':1700,'lipid':1000,'genophore':1560,'virus':1381}#Filipe: 'virus':1455,

# Physical constants [SI-units]
DICT_physical_constants = {'e':1.602176487E-19,'c':299792458,'h':6.62606896E-34,'re':2.8179402894E-15,'barn':1E-28,'u':1.66053886E-27}

# Load scattering factors and atomic masses from file
DICT_atomic_mass = {}
DICT_scattering_factors = {}
unpickle_scattering_factors()

# Direct output to OUT object.
commandline_out_activate()
