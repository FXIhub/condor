#================================#
# Python tools - Multiprocessing #
#================================# 
#
# Author: Max Hantke
# Email: maxhantke@gmail.com

from pylab import *
import multiprocessing
import time

def worker_call(worker,pipe):
    work_package = pipe.recv()
    while work_package != "fika":
        result = worker(work_package)
        pipe.send(result)
        work_package = pipe.recv()        

def multiprocess(Nprocesses,Njobs,worker,getwork,logres=None,logger=None):
    pipes_end_host = list(zeros(Nprocesses))
    pipes_end_worker = list(zeros(Nprocesses))
    processes = list(zeros(Nprocesses))

    for i in range(Nprocesses):
        pipes_end_host[i], pipes_end_worker[i] =  multiprocessing.Pipe()
        processes[i] = multiprocessing.Process(target=worker_call,
                                               args=(worker,pipes_end_worker[i],) )
        processes[i].start()
        if logger != None:
            logger.debug("Process %i started",i)

    Njobs_done = 0
    Njobs_started = 0

    # Send out initial jobs
    for i in range(Nprocesses):
        pipes_end_host[i].send(getwork())
        Njobs_started+=1
        if logger != None:
            logger.debug("Initial job for process %i sent",i)

    t_start = time.time()

    while Njobs_done < Njobs:
        for i in range(Nprocesses):
            if pipes_end_host[i].poll():
                result = pipes_end_host[i].recv()
                if logger != None:
                    logger.info("Datarate %.1f Hz; job %i/%i; process %i/%i",Njobs_done/(time.time()-t_start),Njobs_done,Njobs,i,Nprocesses)
                Njobs_done += 1
                if logres != None:
                    logres(result)
                if Njobs_started < Njobs:
                    pipes_end_host[i].send(getwork())
                    Njobs_started += 1


        time.sleep(0.1)

    for i in range(Nprocesses):
        pipes_end_host[i].send("fika")
        processes[i].join()
        pipes_end_host[i].close()

# for testing
def my_worker(pipe):
    seed()
    A = randint(10)
    B = randint(10)
    res =  "%i+%i=%i" % (A,B,A+B)
    print res
    return res

def my_getwork():
    return randint(10),randint(10)

