/*
 * WtConvol.c
 * A function used in HMMs for molecular motors.
 *
 * Form the product z[w] = sum_u { a[u] * b[u+w] * m[u,u+w] }  
 *  where the following sizes are assumed:
 *        a is nu x 1,
 *        b is nu x 1,
 *        m is nu x nu.
 *
 * This is a MEX-file for MATLAB.
 * Matlab call: z=WeightedConvol(a,b,m)
 *
 * F. Sigworth
 * 4 April 2007
 *
 */

#include "mex.h"

void WtConvol(double *a, double *b, double *m, double *z, mwSize nu)
{
  mwSize u,w;
  double sum;
  double *aptr, *bptr, *mptr;
  
// Simple code that runs about the same speed as the obscure code below.
  for (w=0; w<nu; w++) {
      
      sum=0.0;      
      for (u=0; u<nu-w; u++) {
        sum+=a[u] * b[u+w] * m[u+w+u*nu];
      }
      for (u=nu-w; u<nu; u++) {
        sum+=a[u] * b[u+w-nu] * m[u+w+(u-1)*nu];
      }
    z[w]=sum;
  }
 }

//   for (w=0; w<nu; w++) {
//       
//       //  we accumulate into z[w] the product
//       aptr=a;      //   a[u]
//       bptr=b+w;    // * b[u+w]
//       mptr=m+w;    // * m[u+w,u]
//       sum=0.0;
//       // the sum is accumulated in two parts, to include wrap-around.
//       for (u=0; u<nu-w; u++) {
//         sum+=*(aptr++) * *(bptr++) * *mptr;
//         mptr+=nu+1;  // shift over one column with each row
//       }
//       bptr-=nu;     // wrap to the top row
//       mptr-=nu;     // wrap
//       for (u=nu-w; u<nu; u++) {
//         sum+=*(aptr++) * *(bptr++) * *mptr;
//         mptr+=nu+1;  // shift over one column
//       }
//     *(z++)=sum;
//   }
// }


/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  double *a,*b,*m,*z;
  double  x;
  mwSize asize, bsize, msize;
  
  /*  check for proper number of arguments */
  /* NOTE: You do not need an else statement when using mexErrMsgTxt
     within an if statement, because it will never get to the else
     statement if mexErrMsgTxt is executed. (mexErrMsgTxt breaks you out of
     the MEX-file) */
  if(nrhs!=3) 
    mexErrMsgTxt("Three inputs required.");
  if(nlhs!=1) 
    mexErrMsgTxt("One output required.");
  
  /* check the first input argument */
  if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ) {
    mexErrMsgTxt("Input a must be real double vector.");
  }
  asize=mxGetNumberOfElements(prhs[0]);
  a = mxGetPr(prhs[0]);
  
  /* check the second input argument */
  if( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) ) {
    mexErrMsgTxt("Input b must be real double vector.");
  }
  bsize=mxGetNumberOfElements(prhs[1]);
  if ( bsize != asize ) {
    mexErrMsgTxt("Input b must be same size as a.");
  }
  b = mxGetPr(prhs[1]);

  /* check the third input argument */
  if( !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) ) {
    mexErrMsgTxt("Input m must be real double matrix.");
  }
  msize=mxGetNumberOfElements(prhs[2]);
  if ( msize != asize*asize ) {
    mexErrMsgTxt("Input m must be an (nu x nu) matrix.");
  }
  m = mxGetPr(prhs[2]);

  /*  set the output pointer to the output vector */
  plhs[0] = mxCreateDoubleMatrix(asize,1, mxREAL);
  
  /*  create a C pointer to a copy of the output matrix */
  z = mxGetPr(plhs[0]);
  
  /*  call the C subroutine */
  WtConvol(a,b,m,z,asize);
  
}
