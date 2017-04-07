import condor
import h5writer
from mpi4py import MPI

# MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
is_master = (rank == 0)
nworkers = (size - 1)
assert (nworkers > 0), "Need at least 2 MPI processes"

# CXI writer
W = h5writer.H5WriterMPISW(filename='condor.cxi', comm=comm)

# Condor simulation
N = 100
n = N / nworkers
r = N - ((nworkers) * n)
nframes = nworkers * [n]
nframes[0] += r

if not is_master:
    src = condor.Source(wavelength=.8E-9, pulse_energy=1E-3, focus_diameter=1E-6)
    det = condor.Detector(distance=0.58, pixel_size=.6E-3, nx=256, ny=256)
    par = condor.ParticleSphere(diameter=100E-9, material_type="poliovirus")
    E = condor.Experiment(src, {"particle_sphere" : par}, det)
    for i in range(nframes[rank-1]):
        res = E.propagate()
        W.write_slice(res)
W.close(barrier=True)
