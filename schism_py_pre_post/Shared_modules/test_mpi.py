from mpi4py import MPI

# Initialize MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# Print rank and size
print(f"Hello from rank {rank} of {size} processes!")

# Synchronize all processes
comm.Barrier()

# Finalize MPI
MPI.Finalize()