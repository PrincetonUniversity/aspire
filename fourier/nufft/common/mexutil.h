#ifndef __MEXUTIL_H__

#define __MEXUTIL_H__

/*
 * MATLAB utlity routines
 */

// For compatibility with Windows
#define round(x) (x<0?ceil((x)-0.5):floor((x)+0.5))

bool GetDouble(const mxArray *a, double *value)
{
    if( !mxIsDouble(a) || mxIsSparse(a) || mxIsComplex(a))
        return false;
    if(mxGetM(a)*mxGetN(a) > 1)
        return false;
    double *pr = mxGetPr(a);
    
    *value = pr[0];
    return true;
}


bool GetInteger(const mxArray *a, int *value)
{
    if( !mxIsDouble(a) || mxIsSparse(a) || mxIsComplex(a))
        return false;
    if(mxGetM(a)*mxGetN(a) > 1)
        return false;
    double *pr = mxGetPr(a);
    
    // check to see that the value is actually an integer
    if( floor(pr[0])!=pr[0])
        return false;
    
    *value = (int)pr[0];
    return true;
}

bool IsDoubleArray(const mxArray *arr)
{
    if(!mxIsDouble(arr) | mxIsSparse(arr) | mxIsComplex(arr)) {
        return 0;
    }
    return 1;
    
}

bool IsComplexArray(const mxArray *arr)
{
    if(!mxIsComplex(arr) | mxIsSparse(arr)) {
        return 0;
    }
    return 1;
    
}


#endif

