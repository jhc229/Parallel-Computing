#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>

#ifdef __APPLE__
#include <OpenCL/OpenCl.h>
#else
#include <CL/cl.h>
#endif

/* 
   Compiling on HokieSpeed:
   
   module purge
   module load cuda/6.5.14
   module load gcc/5.1.0

   gcc -I/opt/apps/cuda/6.5.14/include/ -o reduce reduce.c -lOpenCL -lm
   
   To run on HokieSpeed:

   ./reduce

*/


void pfn_notify(const char *errinfo, const void *private_info, size_t cb, void *user_data)
{
  fprintf(stderr, "OpenCL Error (via pfn_notify): %s\n", errinfo);
}

int main(int argc, char **argv){

  int plat = 0;
  int dev  = 0;

  /* set up CL */
  cl_double            err;
  cl_platform_id    platforms[100];
  cl_uint           platforms_n;
  cl_device_id      devices[100];
  cl_uint           devices_n ;

  cl_context        context;
  cl_command_queue  queue;
  cl_device_id      device;

  /* get list of platform IDs (platform == implementation of OpenCL) */
  clGetPlatformIDs(100, platforms, &platforms_n);

  if( plat > platforms_n) {
    printf("ERROR: platform %d unavailable \n", plat);
    exit(-1);
  }

  // find all available device IDs on chosen platform (could restrict to CPU or GPU)
  cl_uint dtype = CL_DEVICE_TYPE_ALL;
  clGetDeviceIDs( platforms[plat], dtype, 100, devices, &devices_n);

  printf("devices_n = %d\n", devices_n);

  if(dev>=devices_n){
    printf("invalid device number for this platform\n");
    exit(0);
  }

  // choose user specified device
  device = devices[dev];

  // make compute context on device
  context = clCreateContext((cl_context_properties *)NULL, 1, &device, &pfn_notify, (void*)NULL, &err);

  // create command queue
  queue   = clCreateCommandQueue(context, device, CL_QUEUE_PROFILING_ENABLE, &err);

  // build kernel function
  const char *sourceFileName = "oclKernels.cl";
  const char *functionName = "oclIterateKernel";
  const char *functionName2 = "oclReductionKernel";

  // read in text from source file

  struct stat statbuf;
  FILE *fh = fopen(sourceFileName, "r");
  if (fh == 0){
    printf("Failed to open: %s\n", sourceFileName);
    exit(-1);
  }
  /* get stats for source file */
  stat(sourceFileName, &statbuf);

  /* read text from source file and add terminator */
  char *source = (char *) malloc(statbuf.st_size + 1);
  fread(source, statbuf.st_size, 1, fh);
  source[statbuf.st_size] = '\0';

  /* create program from source */
  cl_program program = clCreateProgramWithSource(context, 1, (const char **) & source, (size_t*) NULL, &err);

  if (!program){
    printf("Error: Failed to create compute program!\n");
    exit(-1);
  }

  /* compile and build program */
  const char *allFlags = " ";
  err = clBuildProgram(program, 1, &device, allFlags, (void (*)(cl_program, void*))  NULL, NULL);

  /* check for compilation errors */
  char *build_log;
  size_t ret_val_size;
  err = clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, 0, NULL, &ret_val_size);

  build_log = (char*) malloc(ret_val_size+1);
  err = clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, ret_val_size, build_log, (size_t*) NULL);

  /* to be carefully, terminate with \0
     there's no information in the reference whether the string is 0 terminated or not */
  build_log[ret_val_size] = '\0';

  /* print out compilation log */
  fprintf(stderr, "%s", build_log );

  /* create runnable kernel */
  cl_kernel kernel = clCreateKernel(program, functionName, &err);
  cl_kernel kernel2 = clCreateKernel(program, functionName2, &err);
  if (! kernel || err != CL_SUCCESS){
    printf("Error: Failed to create compute oclIterateKernel kernel!\n");
    exit(-1);
  }
  else if (! kernel2 || err != CL_SUCCESS){
    printf("Error: Failed to create compute oclReductionKernel kernel2!\n");
    exit(-2);
  }



  //int N = 2560; /* vector size  */

  /* read N from the command line arguments */
  int N = atoi(argv[1]);
  int M = N; // by default use a square circuit

  /* use for computed epsilon */
  int Ncells = (N+2)*(M+2); // number of cells

  /* create host array */
  size_t sz = Ncells*sizeof(double);
  size_t szEpsilon = ((Ncells+256-1)/256)*sizeof(double)

  double *h_Iold = (double*) malloc(sz);
  double *h_Inew = (double*) malloc(sz);
  double *h_partEpsilon = (double*) malloc(szEpsilon)

  // fill up host array.
  int n;
  for(n=0;n<Ncells;++n){
    h_Iold[n] = 1;
    h_Inew[n] = 0;
  }
  
  // set current
  // int c_idx(int N, int i, int j)
  // idx(i + j*(N+2)) 
  h_Iold[0 + 1*(N+2)] = 1; 
  h_Inew[0 + 1*(N+2)] = 1; 

  h_Iold[i + j*(N+3)] = 1;
  h_Inew[i + j*(N+3)] = 1;

  //double *c_Inew, *c_Iold, *c_partEpsilon;

  /* create device buffer and copy from host buffer */
  cl_mem c_Iold = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sz, h_Iold, &err);
  cl_mem c_Inew = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sz, h_Inew, &err);
 // cl_mem c_partEpsilon = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, szEpsilon, h_partEpsilon &err);

    int dim = 1;
    int Nt = 256;
    int Ng = Nt*((N+Nt-1)/Nt);
    size_t local[3] = {Nt,1,1};
    size_t global[3] = {Ng,1,1};


  /* iterate using the Jacobi method here */
  int it, Nit=40;
  for(it=0;it<Nit;++it){

    /* iterate from Iold to Inew */
    //	cudaIterate(M, N, c_Iold, c_Inew);
    /* now set kernel arguments */
    clSetKernelArg(kernel, 0, sizeof(int), &M);
    clSetKernelArg(kernel, 1, sizeof(int), &N);
    clSetKernelArg(kernel, 2, sizeof(cl_mem), &c_Iold);
    clSetKernelArg(kernel, 3, sizeof(cl_mem), &c_Inew);

    /* CUDA thread dimensions for two-dimensional array of thread-blocks */
    //	dim3 gDim( (N+Titerate-1)/Titerate, (M+Titerate-1)/Titerate);
    //	dim3 bDim( Titerate, Titerate);
    /* invoke CUDA cudaIterateKernel */
    //	cudaIterateKernel <<< gDim, bDim >>> (M, N, c_Iold, c_Inew);
   
    /* queue up kernel */
    clEnqueueNDRangeKernel(queue, kernel, dim, 0, global, local, 0, (cl_event*)NULL, NULL);


    /* iterate from Inew to Iold */
    //	cudaIterate(M, N, c_Inew, c_Iold);
    clSetKernelArg(kernel, 0, sizeof(int), &M);
    clSetKernelArg(kernel, 1, sizeof(int), &N);
    clSetKernelArg(kernel, 2, sizeof(cl_mem), &c_Inew);
    clSetKernelArg(kernel, 3, sizeof(cl_mem), &c_Iold);
    
    /* queue up kernel */
    clEnqueueNDRangeKernel(queue, kernel, dim, 0, global, local, 0, (cl_event*)NULL, NULL);




    /* compute epsilon (change in current) */
    //	epsilon = cudaCalculateEpsilon(M, N, c_Iold, c_Inew, c_partEpsilon, h_partEpsilon);
    /* Cells in circuit */
    int  L = (N+2)*(M+2);

    clSetKernelArg(kernel2, 0, sizeof(int), &L);
    clSetKernelArg(kernel2, 1, sizeof(int), &C_Iold);
    clSetKernelArg(kernel2, 2, sizeof(cl_mem), &c_Inew);
    clSetKernelArg(kernel2, 3, sizeof(cl_mem), &h_partEpsilon);
 
    /* Assume a one-dimensional thread array */
    //	dim3 gDim( (L+Treduction-1)/Treduction); // number of blocks of size Treduction to cover M 
    //	dim3 bDim(Treduction);                   // number of threads per block
    int NgEpsilon = Nt*((L+Nt-2)/Nt);

    /* call kernel to partialReductionKernel */
    //	partialReductionKernel <<< gDim, bDim >>> (L, c_Iold, c_Inew, c_partEpsilon);
    size_t global2[3] = {NgEpsilon,1,1};
    clEnqueueNDRangeKernel(queue, kernel2, dim, 0, global2, local, 0, (cl_event*)NULL, NULL);
    
    /* Copy partially reduced array from DEVICE to HOST */
    //	cudaMemcpy(h_partEpsilon, c_partEpsilon, gDim.x*sizeof(double), cudaMemcpyDeviceToHost); 
    cl_mem c_partEpsilon = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, szEpsilon, h_partEpsilon &err);  

    /* Finish reduction on HOST */
    double epsilon = 0;
    int n;
    for(n=0;n<NgEpsilon;++n){
      epsilon += h_partEpsilon[n];
    }

    epsilon = sqrt(epsilon);

    /* print current residual */
    printf("epsilon = %g\n", epsilon);

  }

  /* blocking read from device to host */
  clFinish(queue);

  /* blocking read to host */
  clEnqueueReadBuffer(queue, c_partEpsilon, CL_TRUE, 0, sz, h_partEpsilon, 0, 0, 0);

  /* print out results */
  for(n=0;n<N/Nt;++n)
    printf("h_partEpsilon[%d] = %g\n", n, h_partEpsilon[n]);

  exit(0);

}