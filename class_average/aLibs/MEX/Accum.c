/*
 * Accum.c
 * A function used in interpolation and gridding
 *
 * Accum(A,B,indices)
 *  Performs the following calculation:
 *        for i=1:numel(A)
 *          B[indices[i]]+=A[i];
 *        end;
 *
 * Note that it is nonstandard in modifying one of the RHS arguments, namely B.
 *
 * This is a MEX-file for MATLAB.
 *
 * F. Sigworth
 * 13 April 2007
 *
 */

#include "mex.h"

void Accum(double *A, double *B, double *Ind, mwSize na, mwSize nb)
// A and ind must be the same size, and ind is bounded by the size of B.
{
    mwSize i,j;
    for (j=0; j<na; j++){
        i=Ind[j]-1;
        if ((i<0)||(i>=nb)){
          mexErrMsgTxt("Array index out of bounds");
        }
        B[i]+=A[j];
    }
}



/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  double *a,*b,*ind;
  mwSize asize, bsize, indsize;
  
  /*  check for proper number of arguments */
  /* NOTE: You do not need an else statement when using mexErrMsgTxt
     within an if statement, because it will never get to the else
     statement if mexErrMsgTxt is executed. (mexErrMsgTxt breaks you out of
     the MEX-file) */
  if(nrhs!=3) 
    mexErrMsgTxt("Three inputs required.");
  
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
  b = mxGetPr(prhs[1]);

  /* check the third input argument */
  if( !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) ) {
    mexErrMsgTxt("Input ind must be real double matrix.");
  }
  indsize=mxGetNumberOfElements(prhs[2]);
  if ( indsize != asize ) {
    mexErrMsgTxt("Input ind must be the same size as A.");
  }
  ind = mxGetPr(prhs[2]);

  /*  call the C subroutine */
  Accum(a,b,ind,asize,bsize);

}
