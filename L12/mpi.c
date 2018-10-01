#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int main(int argc, char **argv){

  int rank,size;
  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  long long int Ninside = 0; // number of random points inside 1/4 circle
  long long int Ntests = 1000000000;
  long long int n;

  long long int *messageOut=(long long int*) malloc(Ntests*sizeof(long long int));
  long long int *messageIn=(long long int*) malloc(Ntests*sizeof(long long int));
  double estpi = 0;

  int loc=size;

  
    srand48(rank);

    for(n=0;n<Ntests;++n){
      double x = drand48();
      double y = drand48();
    
      if(x*x+y*y<1){
        ++Ninside;
      }
    }

    estpi = 4.*(Ninside/(double)Ntests);

    if(rank!=0){

    int messageDest = 0;
    messageOut[0] = Ninside;
      MPI_Send(messageOut,
	       1,
	       MPI_LONG_LONG_INT,
	       messageDest,
	       Ntests,
	       MPI_COMM_WORLD);
    
  }else if(rank==0){

    MPI_Status status;

      for(loc;loc>rank;--loc){
	
	MPI_Recv(messageIn,
		 1,
		 MPI_LONG_LONG_INT,
		 loc,
		 Ntests,
		 MPI_COMM_WORLD,
		 &status);

	MPI_Reduce(&estpi,&messageIn[loc],1,MPI_LONG_LONG_INT,MPI_SUM,0,MPI_COMM_WORLD);
      }
            
  }
  
  if (rank==0) 
    printf("Our estimate of pi is %g \n", estpi);
  
  MPI_Finalize();
  return 0;
}
