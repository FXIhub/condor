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

from __future__ import print_function, absolute_import # Compatibility with python 2 and 3
import numpy

import logging
logger = logging.getLogger(__name__)

from .log import log_and_raise_error,log_warning,log_info,log_debug

class Profile:
    """
    Class for spatial illumination profile

    **Arguments:**
      
      :model (str): See :meth:`condor.utils.profile.Profile.set_model`
    
      :focus_diameter (float): Focus diameter / full-width half maximum in unit meter

    """
    def __init__(self, model, focus_diameter):
        self.set_model(model)
        self.focus_diameter = focus_diameter
        
    def set_model(self,model):
        """
        Set the model

        Args:
    
          :model (str): Model for the spatial illumination profile
        
            *Choose one of the following options:*

              - ``None`` - flat illumination profile of infinite extent and the intensity at every point corresponds to the intensity of a top hat profile (see below)

              - ``\'top_hat\'`` - Circular top hat profile
    
              - ``\'pseudo_lorentzian\'`` - 2D radially symmetrical profile composed of two Gaussian profiles obtained by a fit to a Lorentzian profile (used instead of a real Lorentzian because this function can be normalised)

              - ``\'gaussian\'`` - 2D radially symmetrical Gaussian profile
        """
        if model is None or model in ["top_hat","pseudo_lorentzian","gaussian"]:
            self._model = model
        else:
            log_and_raise_error(logger, "Pulse profile model %s is not implemented. Change your configuration and try again." % model)
            sys.exit(0)

    def get_model(self):
        """
        Return profile model
        """
        return self._model

    def get_radial(self):
        """
        Return radial function of intensity profile
        """
        if self._model is None:
            # we always hit with full power
            p = lambda r: 1. / (numpy.pi * (self.focus_diameter / 2.)**2)
        elif self._model == "top_hat":
            # focus diameter is diameter of circular top hat profile
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

