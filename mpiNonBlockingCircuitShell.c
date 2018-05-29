#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

/* function to convert from (i,j) cell index to linear storage index */
int idx(int N, int i, int j){

  return i + j*(N+2);
}


double allReduce(double data){
  // INPUT: data value from this MPI process to be summed up
  // OUTPUT: sum of the data variables from all MPI processes

//int MPI_Allreduce(const void *sendbuf, void *recvbuf, int count,
//                  MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)

//int MPI_Iallreduce(const void *sendbuf, void *recvbuf, int count,
//                   MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, MPI_Request *request)
  double someValue;
  //MPI_Comm comm = MPI_COMM_WORLD;
 // int rank, size, i = 0;
  //MPI_Comm_rank(MPI_COMM_WORLD, &rank);
 // MPI_Comm_size(MPI_COMM_WORLD, &size);
 // double *someData = (double*) calloc(size, sizeof(double));
  //for (i = 0; i < )

//It is not strictly necessary to compute e at every iteration. Add a command line argument to
//change the number of iterations between computation of.
  MPI_Allreduce(&data, &someValue, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  return someValue;  
 
}


/* function to compute l2 norm (Euclidean norm of difference between Inew and Iold */
double calculateEpsilon(int M, int N, double *Iold, double *Inew){

  double epsilon = 0;

  int i,j;
  for(j=1;j<=M;++j){
    for(i=1;i<=N;++i){
      epsilon += pow(Inew[idx(N,i,j)]-Iold[idx(N,i,j)],2);
    }
  }
  epsilon = allReduce(epsilon);

  // add up this quantity from all processeszdx
  epsilon = sqrt(epsilon);

  return epsilon;
}

/* reference halo exchange from HW02 */
void haloExchange(int M, int N, double *I){

  int rank, size, tag=999;
  MPI_Status status;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  if(rank<size-1)
    MPI_Send(I+idx(N,1,M), N, MPI_DOUBLE, rank+1, tag, MPI_COMM_WORLD);

  if(rank>0)
    MPI_Send(I+idx(N,1,1), N, MPI_DOUBLE, rank-1, tag, MPI_COMM_WORLD);
  
  if(rank>0)
    MPI_Recv(I+idx(N,1,0), N, MPI_DOUBLE, rank-1, tag, MPI_COMM_WORLD, &status);
  
  if(rank<size-1)
    MPI_Recv(I+idx(N,1,M+1), N, MPI_DOUBLE, rank+1, tag, MPI_COMM_WORLD, &status);
}





void startHaloExchange(int M, int N, double *I, MPI_Request *IsendRequests, MPI_Request *IrecvRequests){

  // initiate halo exchange using MPI_Isend and MPI_Irecv 
  int rank, size, tag=999;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Status status;

  // Your Q2a code starts here
  if(rank<size-1)
  	MPI_Isend(I+idx(N,1,M),   N, MPI_DOUBLE, rank+1,   tag, MPI_COMM_WORLD, IsendRequests);
	//MPI_Wait(IsendRequests, &status);

  if(rank>0)
 	MPI_Isend(I+idx(N,1,1),   N, MPI_DOUBLE, rank-1,   tag, MPI_COMM_WORLD, IsendRequests +1);
	//MPI_Wait(IsendRequests, &status);
  if(rank>0)
  	MPI_Irecv(I+idx(N,1,0), N, MPI_DOUBLE, rank-1, tag, MPI_COMM_WORLD, IrecvRequests+1);

  if(rank<size-1)
  	MPI_Irecv(I+idx(N,1,M+1), N, MPI_DOUBLE, rank+1, tag, MPI_COMM_WORLD, IrecvRequests );



  // Your Q2a code ends here
}

void endHaloRecv(MPI_Request *IrecvRequests){

  // wait for halo data recv to complete using MPI_Wait
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // Your Q2b code starts here
  MPI_Status status;

 // MPI_Wait(IrecvRequests, &status);
  // Your Q2b code ends here

  if(rank<size-1)
	MPI_Wait(IrecvRequests , &status);

  if(rank>0)
	MPI_Wait(IrecvRequests +1, &status);
}

void endHaloSend(MPI_Request *IsendRequests){

  // wait for outgoing halo data send to leave the buffer using MPI_Wait
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // Your Q2c code starts here
  MPI_Status status;

  //MPI_Wait(IsendRequests, &status);
  // Your Q2c code ends here

  if(rank<size-1)
	MPI_Wait(IsendRequests , &status);

  if(rank>0)
	MPI_Wait(IsendRequests +1, &status);
}

/* function to update Inew from Iold */
void iterate(int M, int N, double *Iold, double *Inew){

  MPI_Request IsendRequests[2], IrecvRequests[2];
  int i,j;

  // initializes the swap of the top/bottom rows needed by the iterate step
  startHaloExchange(M, N, Iold, IsendRequests, IrecvRequests);

  // process cell updates that do not require halo data from other processes
  for(j=2;j<=M-1;++j)
    for(i=1;i<=N;++i)
      Inew[idx(N,i,j)] = 0.25*(Iold[idx(N,i+1,j)] + Iold[idx(N,i-1,j)] + 
			       Iold[idx(N,i,j+1)] + Iold[idx(N,i,j-1)]);

  // wait for the incoming halo data to arrive
  endHaloRecv(IrecvRequests); 

  // finish update for bottom cells
  j = 1; 
  for(i=1;i<=N;++i)
    Inew[idx(N,i,j)] = 0.25*(Iold[idx(N,i+1,j)] + Iold[idx(N,i-1,j)] + 
			     Iold[idx(N,i,j+1)] + Iold[idx(N,i,j-1)]);

  // finish update for top cells
  j = M;
  for(i=1;i<=N;++i)
    Inew[idx(N,i,j)] = 0.25*(Iold[idx(N,i+1,j)] + Iold[idx(N,i-1,j)] + 
			     Iold[idx(N,i,j+1)] + Iold[idx(N,i,j-1)]);


  // wait for the outgoing halo data buffer to be available for use
  endHaloSend(IsendRequests);

}





/* function to solve for loop currents using Jacobi iterative method */
void solve(int M, int N, double tol){

  /* use for computed epsilon */
  double epsilon;

  double *Inew = (double*) calloc((M+2)*(N+2), sizeof(double));
  double *Iold = (double*) calloc((M+2)*(N+2), sizeof(double));

  /* set batteries based on MPI process rank*/
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // Your Q1 code to set the ghost cells for the two batteries 
  // based on MPI process rank starts here
  if(rank==0){ // bottom left cell 
    Iold[idx(N,0,1)] = 1;
    Inew[idx(N,0,1)] = 1;
  }

  if(rank==size-1){ // top right cell
    printf("cell\n");
    Iold[idx(N,N+1,M)] = 1;
    Inew[idx(N,N+1,M)] = 1;
  }
  // Your Q1 code to set the ghost cells ends here 
  
  
  /* iterate using the Jacobi method here */
  do{

    /* iterate from Iold to Inew */
    iterate(M, N, Iold, Inew);

    /* iterate from Inew to Iold */
    iterate(M, N, Inew, Iold);

    /* compute epsilon (change in current) */
    epsilon = calculateEpsilon(M, N, Iold, Inew);

    /* print current residual */
    if(rank==0)
      printf("epsilon = %g\n", epsilon);

  }while(epsilon>tol);

  /* print out the loop current in cell (1 1) and (10 10) */
  if(rank==0)
    printf("I_{11} = %g\n", Iold[idx(N,1,1)]);
  if(rank*M <= 10 && 10 <= ((rank+1)*M))
    printf("I_{10 10} = %g\n", Iold[idx(N,10, 10-M*rank)]);
  
}

/* usage: ./main 100 1e-6
   mpiexec -n 4 ./main 100 25 1e-6
   to solve for a network of 100x100 to tolerance 1e-6 */

int main(int argc, char **argv){

  // Your Q1 code to call MPI_Init starts here 
  MPI_Init(&argc, &argv);
  // Your Q1 code to call MPI_Init ends here 

  {
    /* read N from the command line arguments */
    int N = atoi(argv[1]);
    
    /* read N from the command line arguments */
    int M = atoi(argv[2]);
    
    /* read the user supplied convergence from the command line arguments */
    double tol = atof(argv[3]);
    
    /* perform Jacobi iteration to solve for loop currents in resistor network */
    solve(M, N, tol);
  }

  // Your Q1 code to call MPI_Finalize starts here 
  MPI_Finalize();
  // Your Q1 code to call MPI_Finalize ends here 

  return 0;

}