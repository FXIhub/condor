import numpy
import logging
logger = logging.getLogger('Condor')

class Variation:
    def __init__(self,mode=None,spread=None,n=None,number_of_dimensions=1,name=""):
        self.name = ""
        self._number_of_dimensions = number_of_dimensions
        self.set_mode(mode,spread,n)
        self.reset_counter()

    def reset_counter(self):
        self._i = 0

    def set_mode(self,mode,spread,n):
        if mode not in [None,"normal","poisson","normal_poisson","uniform","range"]:
            logger.error("Variation object cannot be configured with illegal mode %s",mode)
            return
        if mode in ["normal","normal_poisson","uniform","range"] and spread is None:
            logger.error("Variation object cannot be configured because mode \'%s\' requires valid keyword for \'spread\'",mode)
            return
        if mode in ["range"]:
            if n is None:
                logger.error("Variation object cannot be configured because mode \'%s\' requires valid keyword for \'n\'",mode)
                return
            else:
                if self._number_of_dimensions < 1:
                    logger.error("Variation object does not accept values smaller 1 for \'number_of_dimensions\'.")
                    return
                elif self._number_of_dimensions > 2:
                    logger.error("Variation object does not accept values greater 2 for \'number_of_dimensions\'.")
                    return
                elif self._number_of_dimensions == 1:
                    self._grid = numpy.array([numpy.linspace(-spread/2.,spread/2.,n)])
                elif self._number_of_dimensions == 2:
                    X,Y = numpy.meshgrid(numpy.linspace(-spread[1]/2.,spread[1]/2.,n),numpy.linspace(-spread[0]/2.,spread[0]/2.,n))
                    self._grid = numpy.array([Y.flatten(),X.flatten()])
        else:
            self._grid = None
        if self._number_of_dimensions == 1:
            self._spread = [spread]
        else:
            self._spread = spread
        self._mode = mode
        self._n = n
                
    def get(self,v0):
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
            v1 = numpy.random.normal(numpy.random.poisson(v0),self._spread[dim]/2.)
        elif self._mode == "poisson":
            v1 = numpy.random.poisson(v0)
        elif self._mode == "uniform":
            v1 = numpy.random.uniform(v0-self._spread[dim]/2.,v0+self._spread[dim]/2.)
        elif self._mode == "range":
            v1 = v0 + self._grid[dim,self._i % self._n]
        return v1
