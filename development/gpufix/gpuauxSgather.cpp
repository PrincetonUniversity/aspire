/* Gather a single-precision matrix from the GPU.
 *
 * Yoel Shkolnisky, July 2016.
 */

/* Compile using
 * mex gpuauxSgather.cpp -O -I/usr/local/cuda/targets/x86_64-linux/include/ -L/usr/local/cuda/targets/x86_64-linux/lib/ -lcublas
 */

#include <stdint.h>
#include "mex.h"
#include "cublas.h"

//#define DEBUG

void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    uint64_t gptr; // Pointer to GPU array to return.
    int M,N;       // Dimensions of the returned matrix.
    float* A;      // Elements of the GPU-stored matrix.
    cublasStatus retStatus;
          
    if (nrhs != 3) {
        mexErrMsgTxt("gpuauxfree requires 3 input arguments (gptr,m,n)");
    } else if (nlhs !=1 ) {
        mexErrMsgTxt("gpuauxfree requires 1 output argument (A)");
    }
    
    gptr = (uint64_t) mxGetScalar(prhs[0]);
    M=(int) mxGetScalar(prhs[1]);
    N=(int) mxGetScalar(prhs[2]);
    
    #ifdef DEBUG
    mexPrintf("M=%d  N=%d\n",M,N);
    #endif
    
    plhs[0]=mxCreateNumericMatrix(M, N, mxSINGLE_CLASS, mxREAL);
    A=(float*)mxGetPr(plhs[0]);
    
    retStatus = cublasInit();
    if (retStatus != CUBLAS_STATUS_SUCCESS) {
        mexPrintf("CUBLAS: an error occured in cublasInit\n");
    }
    #ifdef DEBUG
    else {
        mexPrintf("CUBLAS: cublasInit worked\n");
    }
    #endif
    
    retStatus = cublasGetMatrix (M, N, sizeof(float), (void*)gptr, M, A, M);
    if (retStatus != CUBLAS_STATUS_SUCCESS) {
        mexPrintf("CUBLAS: an error occured in cublasGetMatrix\n");
    } 
    #ifdef DEBUG
    else {
        printf("CUBLAS: cublasGetMatrix worked\n");
    }
    #endif
    
    cublasShutdown();
    
    #ifdef DEBUG
    mexPrintf("Matrix retrieved from GPU\n");
    #endif
}
