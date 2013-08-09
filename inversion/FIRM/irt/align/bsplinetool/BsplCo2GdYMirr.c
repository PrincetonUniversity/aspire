/*******************************************************************************

        This file is a mex file to wrap C functions of B-spline
        gradient interpolation using mirror condition.

        Coded by Se Young Chun, Oct 30, 2007, the University of Michigan

*******************************************************************************/

#include "mex.h"
#include "matrix.h"

#include "batchinterpol2D.h"
#include "batchinterpol3D.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int dim, m, n, l, ntot, dims[3];
 	float *data1, *data2, *wx, *wy, *wz;
	double spline, nx, ny, nz, offx, offy, offz, mx, my, mz;
	size_t sizebuf;

	if (nrhs == 6)
	{
		spline = mxGetScalar(prhs[5]);
	}
	else if (nrhs == 5)
	{
		spline = 3;
	}
	else
	{
                mexPrintf("VAL = BsplCo2GdYMirr(COEFF, [nx ny (nz)], [offx "
			"offy (offz)], [mx my (mz)], {wx wy (wz)}, deg);\n");
                mexPrintf("[in]\n");
                mexPrintf("\tCOEFF              : bspline coefficients - "
			"single type\n");
                mexPrintf("\t[nx ny (nz)]       : output dimension\n");
                mexPrintf("\t[offx offy (offz)] : offset of the origin\n");
                mexPrintf("\t[mx my (mz)]       : magnification factor for "
			"each direction\n");
                mexPrintf("\t{wx wy (wz)}       : deformed value vectors\n");
                mexPrintf("\tdeg                : (opt) spline basis degree "
			"{default: 3}\n");
                mexPrintf("[out]\n");
                mexPrintf("\tVAL                : interpolated grad y "
			"values\n");
		return;
	}

    	if (mxIsSingle(prhs[0]) != 1)
        {
                mexErrMsgTxt("First argument must be of single type\n");
                return;
        }

    	/* Retrieve the input data */
    	data1 = mxGetData(prhs[0]);

	if (mxGetNumberOfDimensions(prhs[0]) == 2)
	{
		m = mxGetDimensions(prhs[0])[0];
		n = mxGetDimensions(prhs[0])[1];
		l = 1;
	}
	else if (mxGetNumberOfDimensions(prhs[0]) == 3)
	{
		m = mxGetDimensions(prhs[0])[0];
		n = mxGetDimensions(prhs[0])[1];
		l = mxGetDimensions(prhs[0])[2];
	}
	else 
	{
                mexErrMsgTxt("First argument must be 2D or 3D\n");
		return;
	}

	ntot = mxGetNumberOfElements(prhs[1]);
	dim = ntot;
	sizebuf = mxGetElementSize(prhs[1]);
	if (ntot == 2)
	{
		nx = *((double*)mxGetData(prhs[1]));
		ny = *((double*)(mxGetData(prhs[1]) + sizebuf));
		dims[0] = (int)nx;
		dims[1] = (int)ny;
	}
	else if (ntot == 3)
	{
		nx = *((double*)mxGetData(prhs[1]));
		ny = *((double*)(mxGetData(prhs[1]) + sizebuf));
		nz = *((double*)(mxGetData(prhs[1]) + 2*sizebuf));
		dims[0] = (int)nx;
		dims[1] = (int)ny;
		dims[2] = (int)nz;
	}
	else
	{
                mexErrMsgTxt("Second argument should be either [nx ny nz]"
			" or [nx ny]\n");
		return;
	}

	ntot = mxGetNumberOfElements(prhs[2]);
	sizebuf = mxGetElementSize(prhs[2]);
	if (ntot == 2)
	{
		offx = *((double*)mxGetData(prhs[2]));
		offy = *((double*)(mxGetData(prhs[2]) + sizebuf));
	}
	else if (ntot == 3)
	{
		offx = *((double*)mxGetData(prhs[2]));
		offy = *((double*)(mxGetData(prhs[2]) + sizebuf));
		offz = *((double*)(mxGetData(prhs[2]) + 2*sizebuf));
	}
	else
	{
                mexErrMsgTxt("Third argument should be either [offx offy"
			" offz] or [offx offy]\n");
		return;
	}

	ntot = mxGetNumberOfElements(prhs[3]);
	sizebuf = mxGetElementSize(prhs[3]);
	if (ntot == 2)
	{
		mx = *((double*)mxGetData(prhs[3]));
		my = *((double*)(mxGetData(prhs[3]) + sizebuf));
	}
	else if (ntot == 3)
	{
		mx = *((double*)mxGetData(prhs[3]));
		my = *((double*)(mxGetData(prhs[3]) + sizebuf));
		mz = *((double*)(mxGetData(prhs[3]) + 2*sizebuf));
	}
	else
	{
                mexErrMsgTxt("Fourth argument should be either [mx my mz]"
			" or [mx my]\n");
		return;
	}

	ntot = mxGetNumberOfElements(prhs[4]);
	if (ntot == 0)
	{
	}
	else if (ntot == 2)
	{
		wx = (float*)mxGetData(mxGetCell(prhs[4], 0));
		wy = (float*)mxGetData(mxGetCell(prhs[4], 1));
	}
	else if (ntot == 3)
	{
		wx = (float*)mxGetData(mxGetCell(prhs[4], 0));
		wy = (float*)mxGetData(mxGetCell(prhs[4], 1));
		wz = (float*)mxGetData(mxGetCell(prhs[4], 2));
	}
	else
	{
                mexErrMsgTxt("Fifth argument should be {}, {wx wy wz} or"
			" {wx wy}\n");
		return;
	}

    	/* Create an mxArray for the output data */
	plhs[0] = mxCreateNumericArray(mxGetNumberOfElements(prhs[1]), dims, 
			mxSINGLE_CLASS, mxREAL);

	/* Retrieve the output data */
    	data2 = mxGetData(plhs[0]);

	if (ntot == 0)
	{
		if (dim == 3)
		{
			BatchInterpolatedGradY3DMirr(data2, data1, m, n, l, nx, 
				offx, mx, ny, offy, my, nz, offz, mz, 
				(long)spline); 
		}
		else
		{
			BatchInterpolatedGradY2DMirr(data2, data1, m, n, nx, 
				offx, mx, ny, offy, my, (long)spline); 
		}
	}
	else
	{
		if (dim == 3)
		{
			BatchInterpolatedGradY3DMirrWarp(data2, data1, m, n, l,
 				nx, offx, mx, wx, ny, offy, my, wy, nz, offz, 
				mz, wz, (long)spline); 
		}
		else
		{
			BatchInterpolatedGradY2DMirrWarp(data2, data1, m, n, nx,
 				offx, mx, wx, ny, offy, my, wy, (long)spline); 
		}
	}
}

