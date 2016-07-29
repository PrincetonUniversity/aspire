/* Compute C=A*B using cuda-blas.
 * A, B, and C are matrices of single precision.
 * Yoel Shkolnisky, July 2016.
 * mex complie with:
 * mex gpuauxCmul.cpp -O -I/usr/local/cuda/targets/x86_64-linux/include/ -L/usr/local/cuda/targets/x86_64-linux/lib/ -lcudart -lcublas
 */

#include <stdint.h>
#include <inttypes.h>
#include "mex.h"
//#include "cublas.h"
#include <cuda_runtime.h>
#include "cublas_v2.h"
#include "timings.h"

//#define DEBUG

void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    uint64_t gptr;
    cuComplex *gA,*gB,*gC;
    int mA,nA,mB,nB;
    int transA,transB;
    cuComplex alpha,beta;    
    cublasHandle_t handle;
    cublasStatus_t retStatus;
    cudaError_t cudaStat;

    /* Input:
     *  gA          GPU pointer to A
     *  mA          Number or rows in A
     *  nA          Number of columns in A
     *  gB          GPU pointer to B
     *  mB          Number of rows in B
     *  nB          Number of columns in B
     *  transA      use A or At is the multiplication
     *  transB      use B or Bt is the multiplication
     *
     * Output:
     *  gC          GPU pointer to C     
     */
    
    if (nrhs != 8) {
        mexErrMsgTxt("gpuaxuSmul requires 4 input arguments");
    } else if (nlhs != 1) {
        mexErrMsgTxt("gpuauxSmul requires 1 output argument");
    }
    
    /* Get pointers to matrices on the GPU */
    gptr=(uint64_t) mxGetScalar(prhs[0]);
    gA = (cuComplex*) gptr;
    mA=(int) mxGetScalar(prhs[1]);
    nA=(int) mxGetScalar(prhs[2]);
    
    #ifdef DEBUG
    mexPrintf("[%s,%d] gA=%" PRIu64 " mA=%d  nA=%d\n", __FILE__,__LINE__,(uint64_t)gA, mA, nA);
    #endif
    
    gptr=(uint64_t) mxGetScalar(prhs[3]);
    gB=(cuComplex*) gptr;
    mB=(int) mxGetScalar(prhs[4]);
    nB=(int) mxGetScalar(prhs[5]);

    #ifdef DEBUG
    mexPrintf("[%s,%d] gB=%" PRIu64 " mB=%d  nB=%d\n", __FILE__,__LINE__,(uint64_t)gB, mB, nB);
    #endif
        
    alpha.x = (float)1.0;
    alpha.y = (float)0.0;
    beta.x = (float)0.0;
    beta.y = (float)0.0;

    /* M->mA K->nA  L->mB   N-> nB */
        
    /* STARTUP   CUBLAS */
    TIC;
    retStatus = cublasCreate(&handle);
    if (retStatus != CUBLAS_STATUS_SUCCESS) {
        printf("[%s,%d] an error occured in cublasInit\n",__FILE__,__LINE__);
    } 
    #ifdef DEBUG 
    else {
        printf("[%s,%d] cublasInit worked\n",__FILE__,__LINE__);
    }
    #endif    
    TOCM("init");
    /*
     */
    
    /* Allocate C on the GPU */
    TIC;    
    cudaStat = cudaMalloc((void**)&gC,mA*nB*sizeof(cuComplex));
    TOCM("Allocate C");
    if (cudaStat != cudaSuccess) {
        mexPrintf ("[%s,%d] device memory allocation failed",__FILE__,__LINE__);
        //return EXIT_FAILURE; 
    }
    
    #ifdef DEBUG
    mexPrintf("[%s,%d] gC=%" PRIu64 " mC=%d  nC=%d\n", __FILE__,__LINE__,(uint64_t)gC, mA, nB);
    #endif
        
    TIC;
    /* Multiply */
    (void) cublasCgemm (handle,CUBLAS_OP_N ,CUBLAS_OP_N,mA,nB,nA,&alpha,gA,mA,gB,mB,&beta,gC,mA);
    TOCM("mul");
    #ifdef DEBUG
    mexPrintf("[%s,%d] multiply\n", __FILE__,__LINE__);
    #endif
    
    /* Return pointer */
    plhs[0] = mxCreateNumericMatrix(1, 1, mxUINT64_CLASS, mxREAL);
    (*((uint64_t*)mxGetPr(plhs[0]))) = (uint64_t) gC;
    #ifdef DEBUG
    mexPrintf("[%s,%d] return result\n", __FILE__,__LINE__);
    #endif
    
    
    /* Shutdown */
    TIC;
    cublasDestroy(handle);
    TOCM("shutdown");
    #ifdef DEBUG
    mexPrintf("[%s,%d] shutdown\n", __FILE__,__LINE__);
    #endif

}

