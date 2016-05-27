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

import numpy, time, sys, inspect

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

