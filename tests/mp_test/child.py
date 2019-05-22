import platform
import time
from mpi4py import MPI

# Get child process information
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
machine = platform.node()
x = 0
while x < 10:
    print("In process %s x is %s" % (rank, x))
    time.sleep(3)
    x += 1
