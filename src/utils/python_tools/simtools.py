#===========================#
# Python tools - simulation #
#===========================# 
#
# Author: Max Hantke
# Email: maxhantke@gmail.com

import pylab,math

def random_euler_angles():
    r1,r2,r3 = pylab.random(3)
    q1 = pylab.sqrt(1.0-r1)*pylab.sin(2.0*pylab.pi*r2)
    q2 = pylab.sqrt(1.0-r1)*pylab.cos(2.0*pylab.pi*r2)
    q3 = pylab.sqrt(r1)*pylab.sin(2.0*pylab.pi*r3)
    q4 = pylab.sqrt(r1)*pylab.cos(2.0*pylab.pi*r3)
    phi = math.atan2(2.0*(q1*q2+q3*q4), 1.0-2.0*(q2**2+q3**2))
    theta = math.asin(2.0*(q1*q3-q4*q2))
    psi = math.atan2(2.0*(q1*q4+q2*q3), 1.0-2.0*(q3**2+q4**2))
    return [phi,theta,psi]
