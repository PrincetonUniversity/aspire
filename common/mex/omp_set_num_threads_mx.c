#include <mex.h>

#include <math.h>
#include <omp.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int n_threads;

    if(nrhs != 1)
    {
        mexErrMsgTxt("Incorrect number of arguments.");
    }

    if(!mxIsNumeric(prhs[0]) || mxGetNumberOfElements(prhs[0]) != 1)
    {
        mexErrMsgTxt("n_threads muts be a numeric scalar.");
    }

    if(fabs(mxGetScalar(prhs[0])-round(mxGetScalar(prhs[0]))) > 1e-15 ||
        mxGetScalar(prhs[0]) < 1)
    {
        mexErrMsgTxt("n_threads must be a positive integer.");
    }

    n_threads = round(mxGetScalar(prhs[0]));

    omp_set_num_threads(n_threads);
}
