/*
 * QuickShiftMex.c
 *
 * Perform the circshift function for 2D arrays
 *
 * This is a MEX-file for MATLAB.
 * Matlab call: ms=QuickShift(m,s)
 *
 * F. Sigworth
 * 4 April 2007
 *
 */

#include "mex.h"

void QuickShiftMex(double *m, mwSize *dims, mwSize si, mwSize sj, double *z)
{
   mwSize i,j;
  
  for (j=0; j<dims[1]-sj; j++) { // Loop over the first columns
      
      for (i=0; i<dims[0]-si; i++) {  // Loop over the first rows
        z[i+si+dims[0]*(j+sj)]=m[i+dims[0]*j];
      }
      for (i=dims[0]-si; i<dims[0]; i++) {
        z[i+si+dims[0]*(j+sj-1)]=m[i+dims[0]*j];
      }
  }
  for (j=dims[1]-sj; j<dims[1]; j++) {  // Loop over the wrapped columns
      
      for (i=0; i<dims[0]-si; i++) {
        z[i+si+dims[0]*(j+sj-dims[1])]=m[i+dims[0]*j];
      }
      for (i=dims[0]-si; i<dims[0]; i++) {
        z[i+si+dims[0]*(j+sj-dims[1]-1)]=m[i+dims[0]*j];
      }
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
  double *s,*m,*z;
  mwSize ndims, ssize, si, sj;
  mwSize *dimsizes;
  
  /*  check for proper number of arguments */
  /* NOTE: You do not need an else statement when using mexErrMsgTxt
     within an if statement, because it will never get to the else
     statement if mexErrMsgTxt is executed. (mexErrMsgTxt breaks you out of
     the MEX-file) */
  if(nrhs!=2) 
    mexErrMsgTxt("Two inputs required.");
  if(nlhs!=1) 
    mexErrMsgTxt("One output required.");
  
  /* check the first input argument m */
  if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ) {
    mexErrMsgTxt("Input a must be real double vector.");
  }
  ndims=mxGetNumberOfDimensions(prhs[0]);
  if ( ndims != 2 ) {
      mxErrMsgTxt("Input m must be 2D.");
  }
  dimsizes=mxGetDimensions(prhs[0]);
  m = mxGetPr(prhs[0]);
  
  
  /* check the second input argument shifts */
  if( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) ) {
    mexErrMsgTxt("Input shift must be real double vector.");
  }
  ssize=mxGetNumberOfElements(prhs[1]);
  if ( ssize != 2 ) {
    mexErrMsgTxt("Input s must be a two-element vector.");
  }
  s = mxGetPr(prhs[1]);
  si=s[0];  // Convert to integer
  sj=s[1];

  /*  set the output pointer to the output matrix */
  plhs[0] = mxCreateDoubleMatrix(dimsizes[0],dimsizes[1], mxREAL);
  
  /*  create a C pointer to a copy of the output matrix */
  z = mxGetPr(plhs[0]);
  
  /*  call the C subroutine */
  QuickShiftMex(m,dimsizes,si,sj,z);
  
}
