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
from linalg import rotation

# scattering amplitude from homogeneous spheroid:
# -----------------------------------------------
# Sources:
# Feigin 1987
# Hamzeh,Bragg 1974
#
# a: radius perpendicular to the rotation axis of the ellipsoid (a)
# c: radius along the rotation axis of the ellipsoid (c)
# theta: rotation around x-axis (1st)
# phi: rotation around z-axis (2nd)
#
# F = sqrt(I_0) rho_e p/D r_0 4/3 pi a^2 c [ 3 { sin(qH) - qH cos(qH) } / (qH)^3 ]
#   = sqrt(I_0) rho_e p/D r_0 V f(a,c,theta,phi,qx,qy)
# f = 3 { sin(qH) - qH cos(qH) } / (qH)^3
# H = sqrt(asq sin^2(g)+csq cos^2(g))
# g = arccos( ( -qX cos(theta) sin(phi) + qY cos(theta) cos(phi) ) / sqrt(qX^2+qY^2) )
# K = I_0 (rho_e p/D r_0 V)^2
# S = I_0 rho_e^2 = K / (p/D r_0 V)^2
# ============================================================================================
# I = K [ f(asq,csq,theta,phi,qx,qy) ]^2
# ============================================================================================
_q_spheroid_diffraction = lambda qX,qY: numpy.sqrt(qX**2+qY**2)
_g_spheroid_diffraction = lambda qX,qY,theta,phi: numpy.arccos((-qX*numpy.cos(theta)*numpy.sin(phi)+qY*numpy.cos(theta)*numpy.cos(phi))/(_q_spheroid_diffraction(qX,qY)+numpy.finfo("float64").eps))
_H_spheroid_diffraction = lambda qX,qY,a,c,theta,phi: numpy.sqrt(a**2*numpy.sin(_g_spheroid_diffraction(qX,qY,theta,phi))**2+c**2*numpy.cos(_g_spheroid_diffraction(qX,qY,theta,phi))**2)
_qH_spheroid_diffraction = lambda qX,qY,a,c,theta,phi: _q_spheroid_diffraction(qX,qY)*_H_spheroid_diffraction(qX,qY,a,c,theta,phi)
_F_spheroid_diffraction = lambda K,qX,qY,a,c,theta,phi: numpy.sqrt(abs(K))*3.*(numpy.sin(_qH_spheroid_diffraction(qX,qY,a,c,theta,phi))-_qH_spheroid_diffraction(qX,qY,a,c,theta,phi)*numpy.cos(_qH_spheroid_diffraction(qX,qY,a,c,theta,phi)))/(_qH_spheroid_diffraction(qX,qY,a,c,theta,phi)**3+numpy.finfo("float64").eps)
F_spheroid_diffraction = lambda K,qX,qY,a,c,theta,phi: (_qH_spheroid_diffraction(qX,qY,a,c,theta,phi)**6 < numpy.finfo("float64").resolution)*numpy.sqrt(abs(K)) + (_qH_spheroid_diffraction(qX,qY,a,c,theta,phi)**6 >= numpy.finfo("float64").resolution)*_F_spheroid_diffraction(K,qX,qY,a,c,theta,phi)
_I_spheroid_diffraction = lambda K,qX,qY,a,c,theta,phi: abs(K)*(3.*(numpy.sin(_qH_spheroid_diffraction(qX,qY,a,c,theta,phi))-_qH_spheroid_diffraction(qX,qY,a,c,theta,phi)*numpy.cos(_qH_spheroid_diffraction(qX,qY,a,c,theta,phi)))/(_qH_spheroid_diffraction(qX,qY,a,c,theta,phi)**3+numpy.finfo("float64").eps))**2
I_spheroid_diffraction = lambda K,qX,qY,a,c,theta,phi: (_qH_spheroid_diffraction(qX,qY,a,c,theta,phi)**6 < numpy.finfo("float64").resolution)*abs(K) + (_qH_spheroid_diffraction(qX,qY,a,c,theta,phi)**6 >= numpy.finfo("float64").resolution)*_I_spheroid_diffraction(K,qX,qY,a,c,theta,phi)

def get_spheroid_diffraction_formula(p,D,wavelength,X=None,Y=None):
    if X != None and Y != None:
        qmap = generate_qmap(X,Y,p,D,wavelength)
        qX = qmap[:,:,2]
        qY = qmap[:,:,1]
        qZ = qmap[:,:,0]
        I = lambda K,a,c,theta,phi: I_spheroid_diffraction(K,qX,qY,a,c,theta,phi)
    else:
        qmap = lambda X,Y: generate_qmap(X,Y,p,D,wavelength)
        I = lambda X,Y,K,a,c,theta,phi: I_sphere_diffraction(K,qmap(X,Y)[0],qmap(X,Y)[1],a,c,theta,phi)
    return I

to_spheroid_semi_diameter_a = lambda diameter,flattening: flattening**(1/3.)*diameter/2.
to_spheroid_semi_diameter_c = lambda diameter,flattening: flattening**(-2/3.)*diameter/2.
to_spheroid_diameter = lambda a,c: 2*(a**2*c)**(1/3.)
to_spheroid_flattening = lambda a,c: a/c

def to_spheroid_theta(euler_angle_0,euler_angle_1,euler_angle_2):
    v_z = numpy.array([1.0,0.0,0.0])
    v_y = numpy.array([0.0,1.0,0.0])
    v_rot = rotation(v_y,euler_angle_0,euler_angle_1,euler_angle_2)
    theta = numpy.arcsin(numpy.dot(v_rot,v_z))
    return theta

def to_spheroid_phi(euler_angle_0,euler_angle_1,euler_angle_2):
    v_y = numpy.array([0.0,1.0,0.0])
    v_rot = rotation(v_y,euler_angle_0,euler_angle_1,euler_angle_2)
    v_rot[0] = 0.0
    v_rot = v_rot / numpy.sqrt(v_rot[0]**2+v_rot[1]**2+v_rot[2]**2)       
    phi = numpy.arccos(numpy.dot(v_rot,v_y))
    return phi

def mask_fringe_spheroid(qX,qY,a,c,theta,phi,i_fringe):
    s_mins = get_sphere_diffraction_extrema_positions(i_fringe+1)[0]
    H = _H_spheroid_diffraction(qX,qY,a,c,theta,phi)
    q = numpy.sqrt(qX**2+qY**2)
    s = q*H/2./numpy.pi
    M = s < s_mins[i_fringe]
    if i_fringe != 0:
        M *= s > s_mins[i_fringe-1]
    return M

def mask_min_spheroid(qX,qY,a,c,theta,phi,i_min,ds=0.25):
    s_mins = get_sphere_diffraction_extrema_positions(i_min+1)[0]
    H = _H_spheroid_diffraction(qX,qY,a,c,theta,phi)
    q = numpy.sqrt(qX**2+qY**2)
    s = q*H/2./numpy.pi
    M = (s > s_mins[i_min]-ds)*(s < s_mins[i_min]+ds)
    return M

def spheroid_support(N,dx,a,c,phi):
    Y,X = numpy.indices((N,N))
    X = numpy.float64(X-N/2)*dx
    Y = numpy.float64(Y-N/2)*dx
    T = numpy.arctan2(Y,X)-phi
    R = numpy.sqrt(X**2+Y**2)
    Rmax = a*c/numpy.sqrt((a*numpy.sin(T))**2+(c*numpy.cos(T))**2)
    M = R<=Rmax
    return M               
