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

from .scattering_vector import generate_qmap

_F_sphere_diffraction = lambda K,q,r: numpy.sqrt(abs(K))*3*(numpy.sin(q*r)-q*r*numpy.cos(q*r))/((q*r)**3+numpy.finfo("float64").eps)
F_sphere_diffraction = lambda K,q,r: ((q*r)**6 < numpy.finfo("float64").resolution)*numpy.sqrt(abs(K)) + ((q*r)**6 >= numpy.finfo("float64").resolution)*_F_sphere_diffraction(K,q,r)
r"""
Scattering amplitude from homogeneous sphere (ref. [Feigin1987]_)

.. math::

  F(q) = \sqrt{K} \cdot f(q)

  f(q) =  \frac{ 3 \left[ \sin(qr) - qr \cos(qr) \right]}{ (qr)^3 }

:math:`I_0`: Primary intensity on the sample in unit number of photons per square meter

:math:`\rho_e`: Electron density in unit number of electrons per cubic meter

:math:`p`: Pixel size (i.e. edge length) in unit meter

:math:`D`: Detector distance in unit meter

:math:`r_0`: Classical electron radius in unit meter

:math:`V_r`: Sphere volume in unit cubic meter

Args:
  :K (float): Intensity scaling factor :math:`K = I_0 \left(\rho_e \frac{p}{D} r_0 V_r\right)^2`

  :q (float/array): :math:`q`: Length of scattering vector in unit inverse meter

  :r (float): :math:`r`: Sphere radius in unit meter
"""

_I_sphere_diffraction = lambda K,q,r: abs(K)*(3*(numpy.sin(q*r)-q*r*numpy.cos(q*r))/((q*r)**3+numpy.finfo("float64").eps))**2
I_sphere_diffraction = lambda K,q,r: ((q*r)**6 < numpy.finfo("float64").resolution)*abs(K) + ((q*r)**6 >= numpy.finfo("float64").resolution)*_I_sphere_diffraction(K,q,r)
r"""
Scattering Intensity from homogeneous sphere

.. math::

  I(q) = \left|F(q)\right|^2

Args:
  :K (float): See :func:`condor.utils.sphere_diffraction.F_sphere_diffraction`

  :q (float/array): See :func:`condor.utils.sphere_diffraction.F_sphere_diffraction`

  :r (float): :math:`r`: See :func:`condor.utils.sphere_diffraction.F_sphere_diffraction`
"""

#Fringe_sphere_diffraction = None

#def get_sphere_diffraction_formula(p,D,wavelength,X=None,Y=None):
#    if X != None and Y != None:
#        q = generate_absqmap(X,Y,p,D,wavelength)
#        I = lambda K,r: I_sphere_diffraction(K,q,r)
#    else:
#        q = lambda X,Y: generate_absqmap(X,Y,p,D,wavelength)
#        I = lambda X,Y,K,r: I_sphere_diffraction(K,q(X,Y),r)
#    return I

#def get_sphere_diffraction_extrema_positions(N=1,xtol=1E-12,plot=False):
#    from scipy.optimize import fmin
#
#    f = lambda s:  numpy.sin(2*numpy.pi*s) - 2*numpy.pi*s*numpy.cos(2*numpy.pi*s)
#    f_err = lambda s: f(s)**2
#    g = lambda s:  2*numpy.pi*numpy.cos(2*numpy.pi*s) - 2*numpy.pi*numpy.cos(2*numpy.pi*s) + 4*numpy.pi**2*s*numpy.sin(2*numpy.pi*s)
#    g_err = lambda s: g(s)**2
#
#    s_mins = []
#    s_maxs = []
#    s_guess = 0.7
#    while len(s_mins) < N:
#        s_min = fmin(f_err,s_guess,(),xtol, ftol=0.0001, maxiter=None, maxfun=None, full_output=0, disp=0)
#        s_max = fmin(g_err,s_min+0.2,(),xtol, ftol=0.0001, maxiter=None, maxfun=None, full_output=0, disp=0)
#        s_guess = s_max+0.2
#        s_mins.append(s_min)
#        s_maxs.append(s_max)
#
#    s_mins = numpy.array(s_mins)
#    s_maxs = numpy.array(s_maxs)
#
#    if plot:
#        import matplotlib
#        x = numpy.arange(0,s_mins.max()+1.,0.01)
#        matplotlib.plot(x,f(x))
#        matplotlib.plot(s_mins,numpy.zeros(len(s_mins)),".")
#        matplotlib.plot(s_maxs,numpy.zeros(len(s_maxs)),".")
#        matplotlib.show()
#
#    return [s_mins,s_maxs]

#def mask_fringe_sphere(q,r,i_fringe):
#    s_mins = get_sphere_diffraction_extrema_positions(i_fringe+1)[0]
#    s = q*r/2./numpy.pi
#    M = s < s_mins[i_fringe]
#    if i_fringe != 0:
#        M *= s > s_mins[i_fringe-1]
#    return M

#def mask_min_sphere(q,r,i_min,ds=0.25):
#    s_mins = get_sphere_diffraction_extrema_positions(i_min+1)[0]
#    s = q*r/2./numpy.pi
#    M = (s > s_mins[i_min]-ds)*(s < s_mins[i_min]+ds)
#    return M

