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

import numpy

import logging
logger = logging.getLogger(__name__)

from log import log_and_raise_error,log_warning,log_info,log_debug

class Profile:
    """
    Class for spatial illumination profile
    """
    def __init__(self, model, focus_diameter):
        """
        Initialisation of a Profile instance

          *Choose one of the following options:*

            - ``\'top_hat\'`` - Circular top hat profile

            - ``\'pseudo_lorentzian\'`` - 2D radially symmetrical profile composed of two Gaussian profiles (possible to be normalised) that imitate a Lorentzian profile (``focus_diameter`` represents the full-width half maximum of the Pseudo-Lorentzian profile)

            - ``\'gaussian\'`` - 2D radially symmetrical Gaussian profile (``focus_diameter`` represents the full-width half maximum of the Gaussian profile)

            - ``None`` - infinite extent of the illuminatrion and the intensity at every point corresponds to the intensity within a top hat profile (see above)


        Args:
           :model(str): Model for the spatial illumination profile - either \"top_hat\", \"pseudo_lorentzian\", \"gaussian\" or None (infinite extent of the illuminatrion, intensity same as in case of \"top_hat\")
           :focus_diameter(float): Focus diameter (or characteristic dimension) [m]
        """
        self.set_model(model)
        self.focus_diameter = focus_diameter
        
    def set_model(self,model):
        if model is None or model in ["top_hat","pseudo_lorentzian","gaussian"]:
            self._model = model
        else:
            log_and_raise_error(logger, "Pulse profile model %s is not implemented. Change your configuration and try again.")
            sys.exit(0)

    def get_model(self):
        return self._model

    def get_radial(self):
        if self._model is None:
            # we always hit with full power
            p = lambda r: 1. / (numpy.pi * (self.focus_diameter / 2.)**2)
        elif self._model == "top_hat":
            # focus diameter is diameter of top hat profile
            def p(r):
                if numpy.isscalar(r):
                    return (1.  / (numpy.pi * (self.focus_diameter / 2.)**2)) if r < (self.focus_diameter / 2.) else 0.
                else:
                    return (1.  / (numpy.pi * (self.focus_diameter / 2.)**2)) * (r < (self.focus_diameter / 2.))
        elif self._model == "pseudo_lorentzian":
            # focus diameter is FWHM of lorentzian
            sigma = self.focus_diameter / 2.
            p = lambda r: _pseudo_lorentzian_2dnorm(r, sigma)
        elif self._model == "gaussian":
            # focus diameter is FWHM of gaussian
            sigma = self.focus_diameter / (2.*numpy.sqrt(2.*numpy.log(2.)))
            p = lambda r: _gaussian_2dnorm(r, sigma)
        return p
            
_gaussian = lambda x, sigma: numpy.exp(-x**2/(2*sigma**2))

_gaussian_2dnorm = lambda x, sigma: _gaussian(x, sigma) / ( 2 * numpy.pi * sigma**2 )

_lorentzian = lambda x, sigma: sigma**2 / (x**2 + sigma**2)

_pseudo_lorenzian_A1 = 0.74447313315648778 
_pseudo_lorenzian_A2 = 0.22788162774723308
_pseudo_lorenzian_s1 = 0.73985516665883544
_pseudo_lorenzian_s2 = 2.5588165723260907
_pseudo_lorentzian = lambda x, sigma: _pseudo_lorenzian_A1 * _gaussian(x, _pseudo_lorenzian_s1*sigma) + \
                                      _pseudo_lorenzian_A2 * _gaussian(x, _pseudo_lorenzian_s2*sigma)

_pseudo_lorentzian_2dnorm = lambda x, sigma: _pseudo_lorentzian(x, sigma) / ( 2. * numpy.pi * ( _pseudo_lorenzian_A1 * (_pseudo_lorenzian_s1*sigma)**2 + \
                                                                                                _pseudo_lorenzian_A2 * (_pseudo_lorenzian_s2*sigma)**2 ) )

