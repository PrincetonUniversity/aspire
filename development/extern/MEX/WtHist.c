/*
 * WtHist.c
 * This is a MEX-file for MATLAB.
 * "Weighted Histogram" function for image transformations etc.
 * The mex function is,
 *   function [hist,norm]=WtHist(index,weight,nbins)
 * i.e.
 *   [double[] hist, int32[] norm]=WtHist(int32[] index, double[] weight,
 *   int32 nbins)
 *
 * This implements the equivalent of the following Matlab code:
 *
 * if numel(weight) < numel(index)
 *   error('weight vector smaller than index vector');
 * end;
 * hist=zeros(nbins,1);
 * for i=1:numel(index)
 *  j=index(i);
 *  if (j>0 && j<=nbins)
 *    hist(j)=hist(j)+weight(i);
 *    norm(j)=norm(j)+1;
 *  end;
 * end;
 *
 * The elements of index are checked to be in bounds.  Out-of-bounds
 * elements are not counted.  An error message is issued if the weights
 * vector is smaller than the index vector.
 * This MEX function should be called from
 * a function that does type conversion.
 *
 * F. Sigworth 9 Sep 09
 *
 */

#include "mex.h"


void DoWtHist( long *index,
        double* weight,
        mwSize  numels,
        unsigned long  nbins,
        double* hist,
        long*  norm)
        
{
    long   i;
    unsigned long j;
    
    for (i=0; i<numels; i++) {
        j=index[i]-1;
        if (j<nbins){            
            hist[j]+=weight[i];
            norm[j]+=1;
        }
    }
}
/* the gateway function
 * implements the following:
 * [float[] hist, int32[] norm]
 *         =WtHist(int32[] index, double[] weight, int32 nbins)
 */

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    long*  index;
    mwSize*  nbinp;
    mwSize   nbins;
    mwSize   numels;
    double*  weight;
    double*  hist;
    long*   norm;
    
    /* NOTE: You do not need an else statement when using mexErrMsgTxt
     * within an if statement, because it will never get to the else
     * statement if mexErrMsgTxt is executed. (mexErrMsgTxt breaks you out of
     * the MEX-file) */
    
    if(nrhs != 3)
        mexErrMsgTxt("Needs 3 input arguments.");
    
    if( !mxIsInt32(prhs[0]))
        mexErrMsgTxt("Arg 1 must be int32 array.");
    index=(long*)mxGetPr(prhs[0]);
    numels=mxGetM(prhs[0]);
    
    if( !mxIsDouble(prhs[1]) )
        mexErrMsgTxt("Arg 2 must be a double array.");
    weight = mxGetPr(prhs[1]);
    if( mxGetM(prhs[1])<numels)
        mexErrMsgTxt("Arg 2 'weight' vector is smaller than the 'index' vector.");    
    if( !mxIsInt32(prhs[2]))
        mexErrMsgTxt("Arg 3 'nbins' must be int32.");
    nbinp=(mwSize*)mxGetPr(prhs[2]);
    nbins=*nbinp;
    
    if(nlhs != 2)
        mexErrMsgTxt("Needs 2 output arguments");
    plhs[0] = mxCreateDoubleMatrix(nbins, 1, mxREAL);
    hist=mxGetPr(plhs[0]);
    
    plhs[1] = mxCreateNumericMatrix(nbins, 1,
            mxINT32_CLASS, mxREAL);
    norm=(long*)mxGetPr(plhs[1]);
    
    DoWtHist( index, weight, numels, nbins, hist, norm);
    
    return;
}
