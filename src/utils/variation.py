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
"""
Tools for applying noise / statistical variation to values
"""

import numpy
import collections

import logging
logger = logging.getLogger(__name__)

from log import log_and_raise_error,log_warning,log_info,log_debug

class Variation:
    """
    Class for statistical variation of a variable

    Args:
      :mode (str): See :meth:`condor.utils.variation.Variation.set_mode`

      :spread (float/array): See :meth:`condor.utils.variation.Variation.set_spread`

    Kwargs:
      :n (int): Number of samples within the specified range (default ``None``)

      :number_of_dimensions (int): Number of dimensions of the variable (default ``1``)    
    """
    
    def __init__(self,mode,spread,n=None,number_of_dimensions=1):
        self.set_number_of_dimensions(number_of_dimensions)
        self.set_mode(mode)
        self.set_spread(spread)
        self.n = n
        self.reset_counter()
        self.validate()
        
    def get_conf(self):
        """
        Get configuration in form of a dictionary. Another identically configured Variation instance can be initialised by:

        .. code-block:: python

          conf = V0.get_conf()                          # V0: already existing Variation instance
          V1 = condor.utils.variation.Variation(**conf) # V1: new Variation instance with the same configuration as V0
        """
        conf = {}
        conf["mode"] = self.get_mode()
        conf["spread"] = self.get_spread()
        conf["n"] = self.n
        conf["number_of_dimensions"] = self._number_of_dimensions
        return conf
        
    def reset_counter(self):
        """
        Set counter back to zero

        This counter is relevant only if ``mode=\'range\'``
        """
        self._i = 0

    def set_number_of_dimensions(self, number_of_dimensions):
        if number_of_dimensions < 1 or number_of_dimensions > 3:
            log_and_raise_error(logger, "Number of dimensions for variation objects can be only either 1, 2 or 3.")
            return
        self._number_of_dimensions = number_of_dimensions
        
    def get_number_of_dimensions(self):
        """
        Return the number of dimensions of the variation variable
        """
        return self._number_of_dimensions

    def set_mode(self,mode):
        r"""
        Set the mode of variation

        Args:
          :mode (str): Mode of variation

            *Choose one of the following options*
           
            ======================= ==================================================================== =====================================================
            ``mode``                Variation model                                                      ``spread`` parameter
            ======================= ==================================================================== =====================================================
            ``None``                No variation                                                         -
            ``'uniform'``           Uniform random distribution                                          Width
            ``'normal'``            Normal random distribution                                           Standard deviation
            ``'poisson'``           Poissonian random distribution                                       -
            ``'normal_poisson'``    Normal on top of Poissonian random dist.                             Std. dev. of normal distribution
            ``'range'``             Equispaced grid around mean center position                          Width
            ======================= ==================================================================== =====================================================
        """
        if mode not in [None,"normal","poisson","normal_poisson","uniform","range"]:
            log_and_raise_error(logger, "Variation object cannot be configured with illegal mode %s" % mode)
            return
        self._mode = mode

    def get_mode(self):
        """
        Return the mode of variation
        """
        return self._mode

    def validate(self):
        mode = self.get_mode()
        spread = self.get_spread()
        number_of_dimensions = self.get_number_of_dimensions()
        if mode in ["normal","normal_poisson","uniform","range"] and spread is None:
            log_and_raise_error(logger, "Variation object configuration is invalid because mode \'%s\' requires \'spread\' to be not None." % mode)
        if mode in ["range"]:
            if self.n is None:
                log_and_raise_error(logger, "Variation object cannot be configured because mode \'%s\' requires that \'n\' is not None." % mode)
        if spread is not None:
            if number_of_dimensions != len(self._spread):
                log_and_raise_error(logger, "Specified number of dimensions (%i) and length of spread array (%i) do not match." % (number_of_dimensions, len(spread)))
                
    def _get_grid(self):
        mode = self.get_mode()
        if mode == "range":
            if self._number_of_dimensions == 1:
                return numpy.array([numpy.linspace(-self._spread/2.,self._spread/2.,n)])
            elif self._number_of_dimensions == 2:
                Y,X = numpy.meshgrid(numpy.linspace(-self._spread[0]/2.,self._spread[0]/2.,n),numpy.linspace(-self._spread[1]/2.,self._spread[1]/2.,n),indexing="ij")
                return numpy.array([Y.flatten(),X.flatten()])
            elif self._number_of_dimensions == 3:
                Z,Y,X = numpy.meshgrid(numpy.linspace(-self._spread[0]/2.,self._spread[0]/2.,n),numpy.linspace(-self._spread[1]/2.,self._spread[1]/2.,n),numpy.linspace(-self._spread[2]/2.,self._spread[2]/2.,n),indexing="ij")
                return numpy.array([Z.flatten(),Y.flatten(),X.flatten()])                
        else:
            return None

    def set_spread(self, spread):
        """
        Set spread of variation (standard deviation or full spread of values, see also :meth:`condor.utils.variation.Variation.set_mode`)

        Args:
          :spread (float): Width of the variation
        """
        if spread is None:
            self._spread = None
        elif isinstance(spread, collections.Iterable):
            self._spread = list(spread)
        else:
            self._spread = [spread]
            
    def get_spread(self):
        """
        Get spread of variation
        """
        if self._spread is None or len(self._spread) > 1:
            return self._spread
        else:
            return self._spread[0]
    
    def get(self, v0):
        """
        Get next value(s)

        Args:
          :v0 (float/int/array): Value(s) without variational deviation
        """
        if self._number_of_dimensions == 1:
            v1 = self._get_values_for_one_dim(v0,0)
        else:
            v1 = []
            for dim in range(self._number_of_dimensions):
                v1.append(self._get_values_for_one_dim(v0[dim],dim))
            v1 = numpy.array(v1)
        self._i += 1        
        return v1
        
    def _get_values_for_one_dim(self,v0,dim):
        if self._mode is None:
            v1 = v0
        elif self._mode == "normal":
            v1 = numpy.random.normal(v0,self._spread[dim])
        elif self._mode == "normal_poisson":
            v1 = numpy.random.normal(numpy.random.poisson(v0),self._spread[dim])
        elif self._mode == "poisson":
            v1 = numpy.random.poisson(v0)
        elif self._mode == "uniform":
            v1 = numpy.random.uniform(v0-self._spread[dim]/2.,v0+self._spread[dim]/2.)
        elif self._mode == "range":
            v1 = v0 + self._get_grid()[dim,self._i % self._n]
        return v1
