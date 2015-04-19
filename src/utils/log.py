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

import h5py, os, numpy, time, sys, inspect

import logging
logger = logging.getLogger("Condor")

def log(logcall,message):
    "Automatically log the current function details."
    # Get the previous frame in the stack, otherwise it would
    # be this function!!!
    func = inspect.currentframe().f_back.f_code
    # Dump the message + the name of this function to the log.
    logcall("%s: %s in %s:%i" % (
        message, 
        func.co_name, 
        func.co_filename, 
        func.co_firstlineno
    ))

class CXIWriter:
    def __init__(self,filename,N):
        self.filename = os.path.expandvars(filename)
        self.f = h5py.File(filename,"w")
        self.N = N
    def write(self,d,prefix="",i=-1):
        for k in d.keys():
            name = prefix+"/"+k
            log(logger.debug, "Writing dataest %s" % name)
            if isinstance(d[k],dict):
                if name not in self.f:
                    self.f.create_group(name)
                self.write(d[k],name,i)
            elif k != "i":
                self.write_to_dataset(name,d[k],d.get("i",i))
    def write_to_dataset(self,name,data,i):
        log(logger.debug, "Write dataset %s of event %i." % (name,i))
        if name not in self.f:
            #print name
            t0 = time.time()
            if numpy.isscalar(data):
                if i == -1:
                    s = [1]
                else:
                    s= [self.N]
                t=numpy.dtype(type(data))
                if t == "S":
                    t = h5py.new_vlen(str)
                axes = "experiment_identifier:value"
            else:
                data = numpy.array(data)
                s = list(data.shape)
                ndims = len(s)
                axes = "experiment_identifier"
                if ndims == 2: axes = axes + ":x"
                elif ndims == 3: axes = axes + ":y:x"
                elif ndims == 4: axes = axes + ":z:y:x"
                if i != -1:
                    s.insert(0,self.N)
                t=data.dtype
            self.f.create_dataset(name,s,t)
            self.f[name].attrs.modify("axes",[axes])
            t1 = time.time()
            log(logger.debug, "Create dataset %s within %.1f sec." % (name,t1-t0))
        if i == -1:
            if numpy.isscalar(data):
                self.f[name][0] = data
            else:
                self.f[name][:] = data[:]
        else:
            if numpy.isscalar(data):
                self.f[name][i] = data
            else:
                self.f[name][i,:] = data[:]
    def close(self):
        self.f.close()

