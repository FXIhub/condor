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
import pickle

def load_atomic_scattering_factors(data_dir):
    with open('%s/sf.dat' % data_dir, 'r') as f:
        atomic_scattering_factors = pickle.load(f)
    return atomic_scattering_factors

def load_atomic_masses(data_dir):
    with open('%s/sw.dat' % data_dir, 'r') as f:
        atomic_masses = pickle.load(f)
    return atomic_masses

def load_atomic_numbers(data_dir):
    with open('%s/z.dat' % data_dir, 'r') as f:
        atomic_numbers = pickle.load(f)
    return atomic_numbers

def load_atomic_compositions():
    # Realative atomic compositions of some relevant material types
    atomic_compositions = {
        'water':       { "H" :     2., "C" :     0., "N" :     0., "O" :     1., "P" :     0., "S" :     0. }, # Water H2O
        'protein':     { "H" :    86., "C" :    52., "N" :    13., "O" :    15., "P" :     0., "S" :     1. }, # Bergh et al. 2008: H86 C52 N13 O15 S
        'dna':         { "H" :    11., "C" :    10., "N" :     4., "O" :     6., "P" :     1., "S" :     0. }, # Bergh et al. 2008: H11 C10 N4 O6 P
        'lipid':       { "H" :    69., "C" :    36., "N" :     0., "O" :     6., "P" :     1., "S" :     0. }, # Bergh et al. 2008: H69 C36 O6 P
        'cell':        { "H" :    23., "C" :     3., "N" :     1., "O" :    10., "P" :     0., "S" :     1. }, # Bergh et al. 2008: H23 C3 N O10 S
        'poliovirus':  { "H" :492388., "C" :332652., "N" : 98245., "O" :131196., "P" :  7501., "S" :  2340. }, # Molla et al. 1991: C332652 H492388 N98245 0131196 P7501 S2340
        'styrene':     { "H" :     8., "C" :     8., "N" :     0., "O" :     0., "P" :     0., "S" :     0. }, # Styrene C8H8
        'sucrose':     { "H" :    22., "C" :    12., "N" :     0., "O" :    11., "P" :     0., "S" :     0. }, # Sucrose C12H22O11
    }
    return atomic_compositions
    
def load_mass_densities():
    # Estimated mass densities of relevant material types
    mass_densities = {
        'water':      995., # at 25 C: O'Neil, M.J. (ed.). The Merck Index - An Encyclopedia of Chemicals, Drugs, and Biologicals. Cambridge, UK: Royal Society of Chemistry, 2013., p. 1868
        'protein':   1350., # Bergh et al. 2008
        'dna':       1700., # Bergh et al. 2008
        'lipid':     1000., # Bergh et al. 2008
        'cell':      1000., # Bergh et al. 2008
        'poliovirus':1340., # Dans et al. 1966
        'styrene':    902., # at 25 C: Haynes, W.M. (ed.). CRC Handbook of Chemistry and Physics. 94th Edition. CRC Press LLC, Boca Raton: FL 2013-2014, p. 3-488
        'sucrose':   1581., # at 17 C: Lide, D.R. (ed.). CRC Handbook of Chemistry and Physics. 79th ed. Boca Raton, FL: CRC Press Inc., 1998-1999., p. 3-172
    }
    return mass_densities
