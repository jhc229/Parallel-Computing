#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <openacc.h>
#include <time.h>


/* 
   to compile: gcc -fopenmp -o ompCircuit ompCircuit.c -lm 
   usage: ./ompCircuit 100 
   
   Note 0: this code only requires one command line argument (N)
   Note 1: the code iterates for 40 iterations independent of the solution and circuit size 
   Note 2: uses default number of OpenMP threads
   Reference: Lec 27 0:00:00 ~ 0:25:00
*/



/* function to convert from (i,j) cell index to linear storage index */
#pragma acc routine seq
int idx(int N, int i, int j){

  return i + j*(N+2);
}

/* function to compute l2 norm (Euclidean norm of difference between Inew and Iold */
//Q1: Replace the OpenMP directives in the calculateEpsilon function with appropriate OpenACC
//directives.
double calculateEpsilon(int M, int N, double * restrict Iold, double * restrict Inew){

  double epsilon = 0;

  int i,j;

  //#pragma omp parallel for private(i) reduction(+:epsilon)
  //#pragma acc kernels
#pragma acc data copy(Inew[0:(N+2)*(M+2)]) 
#pragma acc data copy(Iold[0:(N+2)*(M+2)])
#pragma acc parallel loop reduction(+:epsilon)
  //#pragma acc loop gang
  for(j=1;j<=M;++j){
    //#pragma acc loop vector(128)
#pragma acc loop
    for(i=1;i<=N;++i){
      epsilon += pow(Inew[idx(N,i,j)]-Iold[idx(N,i,j)],2);
    }
  }
  epsilon = sqrt(epsilon);

  return epsilon;
}

/* function to update Inew from Iold */
//Q2: Replace the OpenMP directives in the iterate function with appropriate OpenACC directives.
//Q3: Add an OpenACC copy directive at the start of the for loop in the iterate function.
void iterate(int M, int N, double * restrict Iold, double * restrict Inew){

  int i,j;

  //#pragma omp parallel for private(i) 


#pragma acc data copy(Inew[0:(N+2)*(M+2)]) 
#pragma acc data copy(Iold[0:(N+2)*(M+2)])
#pragma acc parallel loop
  for(j=1;j<=N;++j){
#pragma acc loop 
    for(i=1;i<=N;++i){
      Inew[idx(N,i,j)] 
	= 0.25*(Iold[idx(N,i+1,j)] + Iold[idx(N,i-1,j)] + Iold[idx(N,i,j+1)] + Iold[idx(N,i,j-1)]); // + (i==1 && j==1) + (i==N && j==N));
    }
  }
}

/* function to solve for loop currents using Jacobi iterative method */
void solve(int M, int N){

  /* use for computed epsilon */
  double epsilon;

  double *Inew = (double*) calloc((M+2)*(N+2), sizeof(double));
  double *Iold = (double*) calloc((M+2)*(N+2), sizeof(double));

  // set battery ghost current
  Iold[idx(N,0,1)] = 1;
  Inew[idx(N,0,1)] = 1;

  Iold[idx(N,N+1,M)] = 1;
  Inew[idx(N,N+1,M)] = 1;

  /* iterate using the Jacobi method here */
  int it, Nit = 40;
  for(it=0;it<Nit;++it){

    /* iterate from Iold to Inew */
    iterate(M, N, Iold, Inew);

    /* iterate from Inew to Iold */
    iterate(M, N, Inew, Iold);

    /* compute epsilon (change in current) */
    epsilon = calculateEpsilon(M, N, Iold, Inew);

    /* print current residual */
    printf("epsilon = %g\n", epsilon);

  }

  /* print out the loop current in cell (1 1) and (10 10) */
  printf("I_{11} = %g\n", Iold[idx(N,1,1)]);
  printf("I_{10 10} = %g\n", Iold[idx(N,10,10)]);
  
}

int main(int argc, char **argv){

  /* read N from the command line arguments */
  int N = atoi(argv[1]);
  int M = N; // default to square circuit

  double start = clock();

  /* perform Jacobi iteration to solve for loop currents in resistor network of M rows and N columns */
  solve(M, N);

  double end = clock();
  
  /* compute time take to solve */
  double elapsed = (end-start)/CLOCKS_PER_SEC;
  
  printf("elapsed time = %f\n", elapsed);

  return 0;

}
