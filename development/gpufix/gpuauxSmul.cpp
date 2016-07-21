/* Compute C=A*B using cuda-blas.
 * A, B, and C are matrices of single precision.
 * Yoel Shkolnisky, July 2016.
 */

#include <stdint.h>
#include <inttypes.h>
#include "mex.h"
#include "cublas.h"
#include "timings.h"

//#define DEBUG

void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    uint64_t gptr;
    float *gA,*gB,*gC;
    int mA,nA,mB,nB;
    int transA,transB;
    float alpha,beta;    
    cublasStatus retStatus;

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
    gA = (float*) gptr;
    mA=(int) mxGetScalar(prhs[1]);
    nA=(int) mxGetScalar(prhs[2]);
    
    #ifdef DEBUG
    mexPrintf("gA=%" PRIu64 " mA=%d  nA=%d\n", (uint64_t)gA, mA, nA);
    #endif
    
    gptr=(uint64_t) mxGetScalar(prhs[3]);
    gB=(float*) gptr;
    mB=(int) mxGetScalar(prhs[4]);
    nB=(int) mxGetScalar(prhs[5]);

    #ifdef DEBUG
    mexPrintf("gB=%" PRIu64 " mB=%d  nB=%d\n", (uint64_t)gB, mB, nB);
    #endif
        
    alpha = 1.0;
    beta = 0.0;

    /* M->mA K->nA  L->mB   N-> nB */
        
    /* STARTUP   CUBLAS */
    TIC;
    retStatus = cublasInit();
    if (retStatus != CUBLAS_STATUS_SUCCESS) {
        printf("CUBLAS: an error occured in cublasInit\n");
    } 
    #ifdef DEBUG 
    else {
        printf("CUBLAS: cublasInit worked\n");
    }
    #endif    
    TOCM("init");
    /*
     */
    
    /* Allocate C on the GPU */
    TIC;
    cublasAlloc (mA*nB, sizeof(float), (void**)&gC);
    TOCM("Allocate C");

    #ifdef DEBUG
    mexPrintf("gC=%" PRIu64 " mC=%d  nC=%d\n", (uint64_t)gC, mA, nB);
    #endif
    
    
    TIC;
    /* Multiply */
    (void) cublasSgemm ('n','n',mA,nB,nA,alpha,gA,mA,gB,mB,beta,gC,mA);
    TOCM("mul");
    
    /* Return pointer */
    plhs[0] = mxCreateNumericMatrix(1, 1, mxUINT64_CLASS, mxREAL);
    (*((uint64_t*)mxGetPr(plhs[0]))) = (uint64_t) gC;
    
    
    /* Shutdown */
    TIC;
    cublasShutdown();
    TOCM("shutdown");    
}

