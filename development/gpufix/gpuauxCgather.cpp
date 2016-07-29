/* Gather a complex single-precision matrix from the GPU.
 *
 * Yoel Shkolnisky, July 2016.
 */

/* Compile using
 * mex gpuauxCgather.cpp -O -I/usr/local/cuda/targets/x86_64-linux/include/ -L/usr/local/cuda/targets/x86_64-linux/lib/ -lcublas
 */

#include <stdint.h>
#include "mex.h"
#include "cublas.h"

//#define DEBUG

void cuComplex2floats(const cuComplex *in, float *outr, float* outi, int Ntot)
{
    int i;
    for (i = 0; i < Ntot; i++)
    {
        outr[i]=in[i].x;
        outi[i]=in[i].y;
    }
}


void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    uint64_t gptr; // Pointer to GPU array to return.
    int M,N;       // Dimensions of the returned matrix.
    cuComplex* a;  // Elements of the GPU-stored matrix in cuComplex format.
    float *Ar,*Ai; // Real and imaginary parts of returned mxArray.
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
        
    retStatus = cublasInit();
    if (retStatus != CUBLAS_STATUS_SUCCESS) {
        mexPrintf("[%s,%d] an error occured in cublasInit\n",__FILE__,__LINE__);
    }
    #ifdef DEBUG
    else {
        mexPrintf("[%s,%d] cublasInit worked\n",__FILE__,__LINE__);
    }
    #endif
    
    a  = (cuComplex*) mxMalloc(sizeof(cuComplex)*M*N);
    retStatus = cublasGetMatrix (M, N, sizeof(cuComplex), (void*)gptr, M, a, M);
    if (retStatus != CUBLAS_STATUS_SUCCESS) {
        mexPrintf("[%s,%d] an error occured in cublasGetMatrix\n",__FILE__,__LINE__);
    } 
    #ifdef DEBUG
    else {
        printf("[%s,%d] cublasGetMatrix worked\n",__FILE__,__LINE__);
    }
    #endif
    
    plhs[0]=mxCreateNumericMatrix(M, N, mxSINGLE_CLASS, mxCOMPLEX);
    Ar=(float*)mxGetPr(plhs[0]);
    Ai=(float*)mxGetPi(plhs[0]);
    cuComplex2floats(a,Ar,Ai,M*N);
    
    cublasShutdown();
    mxFree(a);
    
    #ifdef DEBUG
    mexPrintf("[%s,%d] Matrix retrieved from GPU\n",__FILE__,__LINE__);
    #endif
}
