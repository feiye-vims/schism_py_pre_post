# include <mpi.h>

int main(int argc, char const *argv[])
{
  int myRank, nRank, source, dest, tag;
  MPI_Status status;
  
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  MPI_Comm_size(MPI_COMM_WORLD, &nRank);
  
  
  return 0;
}
