/* Compute c=max(A) using cuda-blas.
 * A is of single precision.
 */

#include <stdint.h>
#include <inttypes.h>
#include "mex.h"
#include "cublas.h"
#include "timings.h"

#define DEBUG

void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    uint64_t gptr;
    float *gA;
    int mA,nA;
    int max_val_ind = 8, max_val;
    int temp;
   
    cublasStatus retStatus;

    /* Input:
     *  gA          GPU pointer to A
     *  mA          Number or rows in A
     *  nA          Number of columns in A
     *
     * Output:
     *  max_val_ind The index pf the maximal value        
     */
    
    if (nrhs != 3) {
        mexErrMsgTxt("gpuaxuSmul requires 3 input arguments");
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
    
    
    TIC;
    /* find max */
    printf("sizeof(float) = %d\n", sizeof(float));
    max_val_ind = cublasIsamax(mA*nA,gA, 1);
    printf("Max_Val_ind = %d \n", max_val_ind);
    //printf("temp = %d", temp);
    TOCM("max");
    
    /* Return pointer */
    plhs[0] = mxCreateNumericMatrix(1, 1, mxUINT64_CLASS, mxREAL);
    (*((uint64_t*)mxGetPr(plhs[0]))) = (uint64_t) max_val_ind;
    
    
    /* Shutdown */
    TIC;
    cublasShutdown();
    TOCM("shutdown");    
}

