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

from scattering_vector import generate_qmap

_F_sphere_diffraction = lambda K,q,r: numpy.sqrt(abs(K))*3*(numpy.sin(q*r)-q*r*numpy.cos(q*r))/((q*r)**3+numpy.finfo("float64").eps)
F_sphere_diffraction = lambda K,q,r: ((q*r)**6 < numpy.finfo("float64").resolution)*numpy.sqrt(abs(K)) + ((q*r)**6 >= numpy.finfo("float64").resolution)*_F_sphere_diffraction(K,q,r)
r"""
Scattering amplitude from homogeneous sphere [Feigin1987]_

.. math::

  F = \sqrt{K} f(r,q_x,q_y)

  f = 3 \frac{ \sin(qr) - qr \cos(qr) }{ (qr)^3 }

:math:`I_0`: Primary intensity on the sample in unit number of photons per square meter

:math:`\rho_e`: Electron density in unit number of electrons per cubic meter

:math:`p`: Pixel size (i.e. edge length) in unit meter

:math:`D`: Detector distance in unit meter

:math:`r_0`: Classical electron radius in unit meter

:math:`V`: Sphere volume in unit cubic meter



Args:

  :K (float): Intensity scaling factor `K = I_0 \left(\rho_e \frac{p}{D} r_0 V\right)^2`

  :q (float/array): :math:`q`: Length of scattering vector in unit inverse meter

  :r (float): :math:`r`: Sphere radius in unit meter

.. [Feigin1987]: Feigin 1987
"""

_I_sphere_diffraction = lambda K,q,r: abs(K)*(3*(numpy.sin(q*r)-q*r*numpy.cos(q*r))/((q*r)**3+numpy.finfo("float64").eps))**2
I_sphere_diffraction = lambda K,q,r: ((q*r)**6 < numpy.finfo("float64").resolution)*abs(K) + ((q*r)**6 >= numpy.finfo("float64").resolution)*_I_sphere_diffraction(K,q,r)
r"""
Scattering Intensity from homogeneous sphere

.. math::

  I = \left|F\right|^2

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

