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
//#include "cublas.h"
#include <cuda_runtime.h>
#include "cublas_v2.h"


//#define DEBUG

void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    uint64_t gptr; // Pointer to GPU array to free.
    uint64_t hptr;
    cublasHandle_t handle;
    cublasStatus_t retStatus;
    
    if (nrhs != 2) {
        mexErrMsgTxt("gpuauxfree requires 2 input arguments (handle,gptr)");
    } else if (nlhs > 0 ) {
        mexErrMsgTxt("gpuauxfree requires no output arguments");
    }

    /* Get cublas handle */
    hptr=(uint64_t) mxGetScalar(prhs[0]);
    handle = (cublasHandle_t) hptr;

    gptr = (uint64_t) mxGetScalar(prhs[1]);
    #ifdef DEBUG
    mexPrintf("[%s,%d] starting to free GPU array at address %" PRIu64 "\n", __FILE__,__LINE__,gptr);
    #endif
    
    cudaFree((void*)gptr);
    #ifdef DEBUG
    mexPrintf("[%s,%d] memory freed %" PRIu64 "\n", __FILE__,__LINE__,gptr);
    #endif
}
