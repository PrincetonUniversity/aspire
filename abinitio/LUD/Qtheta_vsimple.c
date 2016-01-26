#include <mex.h>
#include <stdio.h>
#include <math.h>
/*
 * Qtheta.c - 
 *
 * Computational function that update theta and S.
 *
 * This is a MEX-file for MATLAB.
 * 
 * Lanhui Wang
 * Jul 20, 2012
 */
 


void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{ 
  mwSize i,j,k,K;  
  double *phi,*C, *mu, *S, *theta;
  double t;
  const mwSize *dims;
  mwSize     ndim;
  
  /* Check for proper number of arguments. */
  if(nrhs!=3) {
    mexErrMsgTxt("Three inputs required.");
  } else if(nlhs>2) {
    mexErrMsgTxt("Too many output arguments.");
  }
  
  ndim = mxGetNumberOfDimensions(prhs[1]);
  dims = mxGetDimensions(prhs[1]);
  K = dims[2];
  
  /* Create matrix for the return argument. */
  plhs[0] = mxCreateDoubleMatrix(2*K,2*K, mxREAL);
  plhs[1] = mxCreateNumericArray(ndim, dims, mxDOUBLE_CLASS, mxREAL);
  /* Assign pointers to each input and output. */
  phi = mxGetPr(prhs[0]);
  C = mxGetPr(prhs[1]);
  mu = mxGetPr(prhs[2]);
  S = mxGetPr(plhs[0]);
  theta = mxGetPr(plhs[1]);
  t=0;
  /* main routine. k*dim1*dim2+j*dim1+i*/
  for (i=0; i<K; i++) {
    for (j=i+1; j<K; j++) {
      t=0;
      for (k=0; k<2; k++) {
      theta[j*2*K+i*2+k] = C[j*2*K+i*2+k]-*mu*(phi[2*j*2*K+2*i+k]*C[i*2*K+j*2]+phi[(2*j+1)*2*K+2*i+k]*C[i*2*K+j*2+1]);
      t = t + theta[j*2*K+i*2+k]* theta[j*2*K+i*2+k];
      }
      t=sqrt(t);
      for (k=0; k<2; k++) {
          if (t>*mu) {
             theta[j*2*K+i*2+k] = theta[j*2*K+i*2+k]/t;
          }
          else {
              theta[j*2*K+i*2+k] = theta[j*2*K+i*2+k]/(*mu);
          }
      S[2*j*2*K+2*i+k] = theta[j*2*K+i*2+k]*C[i*2*K+j*2];
      S[(2*j+1)*2*K+2*i+k] = theta[j*2*K+i*2+k]*C[i*2*K+j*2+1];
      }
    }
  }
}
