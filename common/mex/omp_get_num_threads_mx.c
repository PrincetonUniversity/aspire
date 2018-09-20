#include <mex.h>

#include <math.h>
#include <omp.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int n_threads;

    #pragma omp parallel default(shared)
    {
        int n_parallel = omp_get_num_threads();
        #pragma omp master
        {
            n_threads = n_parallel;
        }
    }

    plhs[0] = mxCreateDoubleScalar(n_threads);
}
