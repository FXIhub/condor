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

import condor

import logging
logger = logging.getLogger(__name__)

log_and_raise_error = lambda logger, message: log(logger, message, lvl="ERROR", exception=RuntimeError, rollback=2)
log_warning = lambda logger, message: log(logger, message, lvl="WARNING", exception=None, rollback=2)
log_info = lambda logger, message: log(logger, message, lvl="INFO", exception=None, rollback=2)
log_debug = lambda logger, message: log(logger, message, lvl="DEBUG", exception=None, rollback=2)

def log(logger, message, lvl, exception=None, rollback=1):
    logcalls = {"ERROR": logger.error,
                "WARNING": logger.warning,
                "DEBUG": logger.debug,
                "INFO": logger.info}
    if lvl not in logcalls:
        print "%s is an invalid logger level." % lvl
        sys.exit(1)
    logcall = logcalls[lvl]
    #logger_condor = logging.getLogger("condor")
    # This should maybe go into a handler
    if (logger.getEffectiveLevel() >= logging.INFO) or rollback is None:
        # Short output
        msg = "%s" % message
    else:
        # Detailed output only in debug mode
        func = inspect.currentframe()
        for r in range(rollback):
            # Rolling back in the stack, otherwise it would be this function
            func = func.f_back
        code = func.f_code
        msg = "%s\n\t=> in \'%s\' function \'%s\' [%s:%i]" % (message,
                                                              func.f_globals["__name__"],
                                                              code.co_name, 
                                                              code.co_filename, 
                                                              code.co_firstlineno)
    logcall("%s:\t%s" % (lvl,msg))
    if exception is not None:
        raise exception(message)
        
def log_execution_time(logger):
    def st_time(func):
        def st_func(*args, **keyArgs):
            t1 = time.time()
            r = func(*args, **keyArgs)
            t2 = time.time()
            try:
                filename = inspect.getsourcefile(func)
                module = inspect.getmodule(func)
                line = inspect.getsourcelines(func)[1]
                loc = "\'%s\' [%s:%i]" % (func.__name__,
                                         filename,
                                         line)
            except TypeError:
                loc = "\'%s\'" % func.__name__
            msg = "Execution time = %.4f sec\n\t=> in function %s" % (t2 - t1,loc)
            log(logger, msg, "DEBUG", exception=None, rollback=None)
            return r        
        return st_func
    return st_time

        
class CXIWriter:
    def __init__(self,filename,N):
        self.filename = os.path.expandvars(filename)
        self.f = h5py.File(filename,"w")
        self.N = N
    def write(self,d,prefix="",i=-1):
        for k in d.keys():
            name = prefix+"/"+k
            log_debug(logger, "Writing dataest %s" % name)
            if isinstance(d[k],dict):
                if name not in self.f:
                    self.f.create_group(name)
                self.write(d[k],name,i)
            elif k != "i":
                self.write_to_dataset(name,d[k],d.get("i",i))
    def write_to_dataset(self,name,data,i):
        log_debug(logger, "Write dataset %s of event %i." % (name,i))
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
            log_debug(logger, "Create dataset %s within %.6f sec." % (name,t1-t0))
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

