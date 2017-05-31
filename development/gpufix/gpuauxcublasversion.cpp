/* Return cuBlas version
 *
 * Yoel Shkolnisky, July 2016.
 */

#include <stdint.h>
#include <inttypes.h>
#include "mex.h"
#include "cublas.h"



void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{

    int version;
    
    cublasGetVersion(&version);
    mexPrintf("version=%d\n",version);
}