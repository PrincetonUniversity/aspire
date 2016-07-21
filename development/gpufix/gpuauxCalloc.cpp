/* Store a complex single precision matrix on the GPU.
 *
 * Yoel Shkolnisky, July 2016.
 */

/* Compile using
 * mex gpuauxCalloc.cpp -O -I/usr/local/cuda/targets/x86_64-linux/include/ -L/usr/local/cuda/targets/x86_64-linux/lib/ -lcublas
 */

#include <stdint.h>
#include <inttypes.h>
#include "mex.h"
#include "cublas.h"

#define DEBUG


void floats2cuComplex( const float *inr, const float *ini, cuComplex *out, int Ntot)
{
    int i;
    
    for (i = 0; i < Ntot; i++)
    {        
        out[i]={inr[i],ini[i]};
    }
}

void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    int M,N;    // Dimensions of the input matrix.
    float *Ar,*Ai;   // Data elements of the input array (real and imaginary parts).
    cuComplex *a;    // Complex data of A in cuComplex layout.
    cuComplex *gA;   // Pointer to the GPU copy of the data of A.
    uint64_t *gptr;  // Address of the GPU-allocated array.
    cublasStatus retStatus;
    
    if (nrhs != 1) {
        mexErrMsgTxt("gpuauxSalloc requires 1 input arguments (matrix A)");
    } else if (nlhs != 1) {
        mexErrMsgTxt("gpuauxSalloc requires 1 output argument");
    }
    
    Ar = (float*) mxGetPr(prhs[0]);  // Get real data of A.
    Ai = (float*) mxGetPi(prhs[0]);  // Get real data of A.
    M = mxGetM(prhs[0]);   // Get number of rows of A.
    N = mxGetN(prhs[0]);   // Get number of columns of A.
        
    #ifdef DEBUG
    mexPrintf("M=%d  N=%d\n",M,N);
    #endif
    
    
    /* Convert A=(Ar,Ai) into cuComplex format */
    a  = (cuComplex*) mxMalloc(sizeof(cuComplex)*M*N);
    floats2cuComplex(Ar,Ai,a, M*N);

    /* ALLOCATE SPACE ON THE GPU */
    cublasAlloc (M*N, sizeof(cuComplex), (void**)&gA);
    
    // test for error
    retStatus = cublasGetError ();   
    if (retStatus != CUBLAS_STATUS_SUCCESS) {
        mexPrintf("CUBLAS: an error occured in cublasAlloc\n");
    } 
    #ifdef DEBUG
    else {
        mexPrintf("CUBLAS: cublasAlloc worked\n");
    }
    #endif
       
    retStatus = cublasSetMatrix (M, N, sizeof(cuComplex),
            a, M, (void*)gA, M);
    
    if (retStatus != CUBLAS_STATUS_SUCCESS) {
        mexPrintf("CUBLAS: an error occured in cublasSetMatrix\n");
    } 
    #ifdef DEBUG
    else {
        mexPrintf("CUBLAS: cublasSetMatrix worked\n");
    }
    #endif
        
    plhs[0] = mxCreateNumericMatrix(1, 1, mxUINT64_CLASS, mxREAL);
    gptr=(uint64_t*) mxGetPr(plhs[0]);
    *gptr=(uint64_t) gA;
    
    cublasShutdown();
    mxFree(a);
    
    #ifdef DEBUG
    mexPrintf("GPU array allocated at address %" PRIu64 "\n", (uint64_t)gA);
    #endif
}

