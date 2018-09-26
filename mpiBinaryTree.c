#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>

int main(int argc, char **argv){

  int rank, size;

  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int messageN = 1;
  int messageTag = 999;

  int *messageOut = (int*) malloc(messageN*sizeof(int));
  int *messageIn = (int*) malloc(messageN*sizeof(int));

  int alive = size;

  while(alive>1){

    if(rank<alive && rank>=alive/2){
      int messageDest = rank - (alive/2);
      //send
      MPI_Send(messageOut,
	       messageN,
	       MPI_INT,
	       messageDest,
	       messageTag,
	       MPI_COMM_WORLD);
    }

    if(rank<alive && rank<alive/2){
      MPI_Status status;
      int messageSource = rank + (alive/2);
      //receive
      MPI_Receive(messageIn,
		  messageN,
		  MPI_INT,
		  messageSource,
		  messageTag,
		  MPI_COMM_WORLD
		  &status);
    }

    alive = alive/2;
  }

  if (rank==0) printf("ehhhh\n");

  MPI_Finalize();
  return 0;
}