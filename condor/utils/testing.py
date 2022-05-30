from __future__ import print_function, absolute_import # Compatibility with python 2 and 3
import numpy

def same_shape(a0, a1):
    s0 = numpy.array(list(a0.shape))
    s1 = numpy.array(list(a1.shape))
    return numpy.array_equal(s0,s1)
