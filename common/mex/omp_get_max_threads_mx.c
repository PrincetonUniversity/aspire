#include <mex.h>

#include <math.h>
#include <omp.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int n_threads;

    n_threads = omp_get_max_threads();

    plhs[0] = mxCreateDoubleScalar(n_threads);
}
