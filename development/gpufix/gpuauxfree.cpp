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

//#define DEBUG

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

    
    gptr = (uint64_t) mxGetScalar(prhs[0]);
    #ifdef DEBUG
    mexPrintf("[%s,%d] starting to free GPU array at address %" PRIu64 "\n", __FILE__,__LINE__,gptr);
    #endif

    
    retStatus = cublasInit();
    if (retStatus != CUBLAS_STATUS_SUCCESS) {
        mexPrintf("[%s,%d] an error occured in cublasInit\n",__FILE__,__LINE__);
    }
    #ifdef DEBUG
    else {
        mexPrintf("[%s,%d] cublasInit worked\n",__FILE__,__LINE__);
    }
    #endif
    
    cublasFree ((void*)gptr);
    #ifdef DEBUG
    mexPrintf("[%s,%d] memory freed %" PRIu64 "\n", __FILE__,__LINE__,gptr);
    #endif
    
    cublasShutdown();
    
    #ifdef DEBUG
    mexPrintf("[%s,%d] GPU array freed at address %" PRIu64 "\n", __FILE__,__LINE__,gptr);
    #endif
}
