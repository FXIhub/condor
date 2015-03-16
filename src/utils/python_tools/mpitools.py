#====================#
# Python tools - MPI #
#====================# 
#
# Author: Max Hantke
# Email: maxhantke@gmail.com

from pylab import *
import time,h5py,sys
from mpi4py import MPI

def multinode_master(comm,Njobs,getwork,logres=None,logger=None):
    comm_size = comm.Get_size()
    t_start = time.time()

    Njobs_done = 0
    Njobs_started = 0
    results = []
    texec = zeros(shape=(comm_size,2))
    # send initial jobs
    for i in range(1,comm_size):
        comm.send(getwork(),dest=i,tag=1)
        texec[i,0] = time.time()
        Njobs_started += 1

    dummybuffer = array(1, dtype='i')        
    status = MPI.Status()
    Dt = {"wait":[],"send":[],"write":[],"read":[],"receive":[],"work":[]}
    while Njobs_done < Njobs:
        t0 = time.time()
        request = comm.Irecv(dummybuffer,MPI.ANY_SOURCE,0)
        if logger != None:
            logger.debug("Waiting for job to finish.")
            w0 = time.time()
        MPI.Request.Wait(request, status)
        if logger != None:
            w1 = time.time()
            logger.debug("Job finished (waiting time: %.1f sec).",w1-w0)
        i_done = status.source
        texec[i_done,1] = t1 = time.time()
        Njobs_done += 1
        if logger != None:
            logger.info("Datarate %.1f Hz job %i/%i rank %i/%i",Njobs_done/(time.time()-t_start),Njobs_done,Njobs,i_done,comm_size)
        if logger != None:
            logger.debug("Waiting for result to be transferred.")
            r0 = time.time()
        result = comm.recv(source=i_done,tag=1)
        if logger != None:
            r1 = time.time()
            logger.debug("Result transferred (time: %.1f)",r1-r0)
        t2 = time.time()
        if logres != None:
            if logger != None:
                logger.debug("Waiting for result to be saved.")
                s0 = time.time()
            logres(result)
            if logger != None:
                s1 = time.time()
                logger.debug("Result saved (time: %.1f).",s1-s0)
        t3 = time.time()
        Dt["work"].append((texec[i_done,1]-texec[i_done,0])/(1.*(comm_size-1)))
        if Njobs_started < Njobs:
            work = getwork()
            t4 = time.time()
            comm.send(work,dest=i_done,tag=1)
            t5 = texec[i_done,0] = time.time()
            Njobs_started += 1
        else:
            t5 = t4 = time.time()
        Dt["wait"].append(t1-t0)
        Dt["receive"].append(t2-t1)
        Dt["write"].append(t3-t2)
        Dt["read"].append(t4-t3)
        Dt["send"].append(t5-t4)
        S = ""
        for k in Dt.keys():
            S += "%s %f " % (k,Dt[k][-1])
        if logger != None:
            logger.info("Speed: " + S)
        sys.stdout.flush()

    for i in range(1,comm_size):
        comm.send("fika",dest=i,tag=1)

    
    return Dt

def write_times(Dt,outfolder):
    filename = outfolder + "/" + "times.h5"
    f = h5py.File(filename,"w")
    for key in Dt.keys():
        f.create_dataset(key,None,None,Dt[key])
    f.close()


def multinode_slave(comm,worker,logger=None):
    work_package = comm.recv(source=0,tag=1)
    dummydata = array(1, dtype='i')
    while work_package != "fika":
        result = worker(work_package)
        comm.Send([dummydata,MPI.INT],0,0)
        comm.send(result,dest=0,tag=1)
        work_package = comm.recv(source=0,tag=1)


