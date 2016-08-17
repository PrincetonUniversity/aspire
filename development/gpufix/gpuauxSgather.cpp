/* Gather a single-precision matrix from the GPU.
 *
 * Yoel Shkolnisky, July 2016.
 */

/* Compile using
 * mex gpuauxSgather.cpp -O -I/usr/local/cuda/targets/x86_64-linux/include/ -L/usr/local/cuda/targets/x86_64-linux/lib/ -lcublas
 */

#include <stdint.h>
#include "mex.h"
//#include "cublas.h"
#include <cuda_runtime.h>
#include "cublas_v2.h"

//#define DEBUG

void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    uint64_t gptr; // Pointer to GPU array to return.
    int M,N;       // Dimensions of the returned matrix.
    float* A;      // Elements of the GPU-stored matrix.
    uint64_t    hptr;
    cublasHandle_t handle;
    cublasStatus_t retStatus;

          
    if (nrhs != 4) {
        mexErrMsgTxt("gpuauxfree requires 4 input arguments (handle,gptr,m,n)");
    } else if (nlhs !=1 ) {
        mexErrMsgTxt("gpuauxfree requires 1 output argument (A)");
    }
    
    hptr=(uint64_t) mxGetScalar(prhs[0]);
    handle = (cublasHandle_t) hptr;

    gptr = (uint64_t) mxGetScalar(prhs[1]);
    M=(int) mxGetScalar(prhs[2]);
    N=(int) mxGetScalar(prhs[3]);
    
    #ifdef DEBUG
    mexPrintf("M=%d  N=%d\n",M,N);
    #endif
    
    plhs[0]=mxCreateNumericMatrix(M, N, mxSINGLE_CLASS, mxREAL);
    A=(float*)mxGetPr(plhs[0]);
       
    retStatus = cublasGetMatrix (M, N, sizeof(float), (void*)gptr, M, A, M);
    if (retStatus != CUBLAS_STATUS_SUCCESS) {
        mexPrintf("[%s,%d] an error occured in cublasGetMatrix\n",__FILE__,__LINE__);
    } 
    #ifdef DEBUG
    else {
        printf("[%s,%d] cublasGetMatrix worked\n",__FILE__,__LINE__);
    }
    #endif

}
