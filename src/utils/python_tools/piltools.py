#====================#
# Python tools - PIL #
#====================# 
#
# Author: Max Hantke
# Email: maxhantke@gmail.com

import os,re,sys,spimage,Image
from pylab import *

def attach(filenamelist,arrangementmap):
    i = 0
    for x in range(0,arrangementmap.shape[1]):
        for y in range(0,arrangementmap.shape[0]):
            filenumber = int(arrangementmap[y,x])
            if not filenumber == -1:
                I_xy = Image.open(filenamelist[filenumber])
                if i == 0:
                    I = Image.new("RGB",(I_xy.size[0]*arrangementmap.shape[1],I_xy.size[1]*arrangementmap.shape[0]),(256,256,256))
                I.paste(I_xy,(x*I_xy.size[0],y*I_xy.size[1]))
                i+=1
    return I
