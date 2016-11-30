#include <stdio.h>
#include <mpi.h>

int main(int argc, char **argv) {
  int size, rank, i;
  char *message = "OIE!";
  MPI_Status st;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  MPI_Bcast(&message, 1, MPI_CHAR, 0, MPI_COMM_WORLD);

  if(rank == 0){
    for(i = 1; i < size; i ++){
      MPI_Recv(&message, 1, MPI_CHAR, i, MPI_ANY_TAG, MPI_COMM_WORLD, &st);
      printf("Master Received Message: %s from %d\n", message, st.MPI_SOURCE);
    }
  }else{
    MPI_Send(&message, 1, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
  }

  MPI_Finalize();
  return 0;
}
