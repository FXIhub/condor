#!/usr/bin/env python
import argparse
import os
import condor
import condor.utils
from condor.utils.log import log_info
import logging
logger = logging.getLogger("condor")
import time, numpy

#if __name__ == "__main__":
def main():
    parser = argparse.ArgumentParser(description='Condor - simulation of single particle X-ray diffraction patterns')
    parser.add_argument('-v', '--verbose', dest='verbose',  action='store_true', help='verbose mode', default=False)
    parser.add_argument('-d', '--debug', dest='debug',  action='store_true', help='debugging mode (even more output than in verbose mode)', default=False)
    parser.add_argument('-t', '--measure-time', dest='measure_time',  action='store_true', help='Measure execution time', default=False)
    parser.add_argument('-r', '--number-of-repetitions', metavar='number_of_repetitions', type=int,
                        help="number of repetitions (for time measurements)", default=1)
    parser.add_argument('-n', '--number-of-patterns', metavar='number_of_patterns', type=int,
                        help="number of patterns to be simulated", default=1)
    args = parser.parse_args()
    if not os.path.exists("./condor.conf"):
        parser.error("Cannot find configuration file \"condor.conf\" in current directory.")
    if args.verbose:
        logger.setLevel("INFO")
    if args.debug:
        logger.setLevel("DEBUG")

    t_exec = []
    t_write = []    

    t0 = time.time()
    
    for i in range(args.number_of_repetitions):
        
        E = condor.experiment.experiment_from_configfile("./condor.conf")
    
        # FOR BENCHMARKING
        #from pycallgraph import PyCallGraph
        #from pycallgraph.output import GraphvizOutput
        #from pycallgraph import Config
        #from pycallgraph import GlobbingFilter
        #config = Config()
        #config.trace_filter = GlobbingFilter(exclude=[
        #'pycallgraph.*',
        #'numpy.*',
        #])
        #with PyCallGraph(output=GraphvizOutput(),config=config):
        
        W = condor.utils.cxiwriter.CXIWriter("./condor.cxi")
        for i in range(args.number_of_patterns):
            t1 = time.time()
            res = E.propagate()
            t2 = time.time()
            W.write(res)
            t3 = time.time()
            t_exec.append(t2 - t1)
            t_write.append(t3 - t2)
        W.close()

    t4 = time.time()
    
    if args.measure_time:
        t_exec = numpy.array(t_exec)
        t_write = numpy.array(t_write)
        print "TIME MEASUREMENT RESULTS IN SECONDS"
        print "Total time: %.3f" % numpy.round(t4 - t0, 3)
        print "Computation time per image: %.3f" % numpy.round(t_exec.mean(), 3)
        print "( Individual computation times: ", numpy.round(t_exec, 3), " )"
        print "Writing time per image: %.3f" % numpy.round(t_write.mean(), 3)
        print "( Individual writing times: ", numpy.round(t_write, 3), " )"
        
        
