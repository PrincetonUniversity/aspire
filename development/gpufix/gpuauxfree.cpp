/* Free matrix on the GPU.
 *
 * Yoel Shkolnisky, July 2016.
 */

/* Compile using
 * mex gpuauxfree.cpp -O -I/usr/local/cuda/targets/x86_64-linux/include/ -L/usr/local/cuda/targets/x86_64-linux/lib/ -lcublas
 */

#include <stdint.h>
#include <inttypes.h>
#include "mex.h"
#include "cublas.h"

#define DEBUG

void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    uint64_t gptr; // Pointer to GPU array to free.
    cublasStatus retStatus;
    
    if (nrhs != 1) {
        mexErrMsgTxt("gpuauxfree requires one input arguments (gptr)");
    } else if (nlhs > 0 ) {
        mexErrMsgTxt("gpuauxfree requires no output arguments");
    }
            
    retStatus = cublasInit();
    if (retStatus != CUBLAS_STATUS_SUCCESS) {
        mexPrintf("CUBLAS: an error occured in cublasInit\n");
    }
    #ifdef DEBUG
    else {
        mexPrintf("CUBLAS: cublasInit worked\n");
    }
    #endif
    
    gptr = (uint64_t) mxGetScalar(prhs[0]);
    cublasFree ((void*)gptr);
    cublasShutdown();
    
    #ifdef DEBUG
    mexPrintf("GPU array freed at address %" PRIu64 "\n", gptr);
    #endif
}
