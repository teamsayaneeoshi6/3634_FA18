#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>

int main(int argc, char **argv){

  int rank, size;

  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int messageN = 10;
  int messageTag = 999;

  int *messageOut = (int*) malloc(messageN*sizeof(int));
  for (int n = 0; n < messageN; ++n){
    messageOut[n] = n;
  }

  int messageDest = (rank + 1)%15;

  MPI_Send(messageOut,
	   messageN,
	   MPI_INT, //data type from MPI
	   messageDest,
	   messageTag,
	   MPI_COMM_WORLD);
  
  if (rank > 0){

    MPI_Status status;

    int *messageIn = (int*) malloc(messageN*sizeof(int));
    int messageSource = rank - 1;

    MPI_Recv(messageIn,
	     messageN,
	     MPI_INT,
	     messageSource,
	     messageTag,
	     MPI_COMM_WORLD,
	     &status);

    if (messageDest == 0){
      printf("owo\n");
    }

  }

  MPI_Finalize();
  return 0;
}
