// nufftpreparemx.cpp - auxiliry function to optimize the performance of 
// the nonequally-spaced FFT.
//
// Precomputes interpolation tables in case we use the same sampling grid 
// over and over again.
//
// call as nufftpreparemx(x,nu,n,m,b,q);
//
// Returns a structure with the precomuted tables, which are then passed to 
// the function that computes the nonequally-spaced FFT
//
//	Returns a structure containing nu, alpha, exp_k_square, exp_k1_factor, 
//  exp_k2_factor,n,m,b,q.
//
//
// Yoel Shkolnisky, January 2008.


//#define DEBUG

#include <mex.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include "mexutil.h"

#define for if(0);else for


void error_msg(const char* msg)
{
    printf("%s: %s",__FILE__,msg);
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
    /* if there are no input arguments, print out help information */
    if(nrhs!=6) {
        error_msg("Not enough parameters.\n");
        return;
    }
    
    #ifdef DEBUG
    FILE *fid=fopen("debug.txt","w");
    #endif
    
    const double PI=atan(1.0)*4;
    
    const mxArray *x=prhs[0];     // double
    const mxArray *nu=prhs[1];    // double
    
    int n,m,q;
    double b;
    
    GetInteger(prhs[2],&n);
    GetInteger(prhs[3],&m);
    GetDouble(prhs[4],&b);
    GetInteger(prhs[5],&q);
    
    if(!IsDoubleArray(x)) {
        error_msg("x must be a double array");
        return;
    }
        
    if(!IsDoubleArray(nu)) {
        error_msg("nu must be a double array");
        return;
    }
    
    // Get pointers to the data   
    const double *x_pr=mxGetPr(x);
    const double *nu_pr=mxGetPr(nu);
    
    const int len_x=mxGetM(x);                          
    const double c1=4*b*PI;
    const double c2=m*n/(2*PI);
    
	// Create tables for fast interpolation

    
    mxArray *alpha = mxCreateDoubleMatrix(len_x,1,mxREAL);
    mxArray *exp_k_square  = mxCreateDoubleMatrix(q+1,1,mxREAL);
    mxArray *exp_k1_factor = mxCreateDoubleMatrix((q+1)*len_x,1,mxREAL);
    mxArray *exp_k2_factor = mxCreateDoubleMatrix((q+1)*len_x,1,mxREAL);
    
    double *alpha_pr = mxGetPr(alpha);
	double *exp_k_square_pr = mxGetPr(exp_k_square);
	double *exp_k1_factor_pr = mxGetPr(exp_k1_factor);
	double *exp_k2_factor_pr = mxGetPr(exp_k2_factor);
	double tmp,tmp1,tmp2,tmp3;
	const  int ofs=(int)(q/2);

    
    for (int i=0;i<len_x;i++) {

		tmp1=(x_pr[i]*c2)*(x_pr[i]*c2)+(x_pr[i+len_x]*c2)*(x_pr[i+len_x]*c2);
		tmp2=nu_pr[i]*nu_pr[i]+nu_pr[i+len_x]*nu_pr[i+len_x];
		tmp3=-2*x_pr[i]*c2*nu_pr[i]-2*x_pr[i+len_x]*c2*nu_pr[i+len_x];
		tmp=tmp1+tmp2+tmp3;

		alpha_pr[i]=exp(-tmp/(4*b))/c1;

		for (int k=(int)-q/2;k<=q/2;k++) {
			exp_k1_factor_pr[(k+ofs)*len_x+i]=exp(-(-k*(2*x_pr[i]*c2-2*nu_pr[i]))/(4*b));
			exp_k2_factor_pr[(k+ofs)*len_x+i]=exp(-(-k*(2*x_pr[i+len_x]*c2-2*nu_pr[i+len_x]))/(4*b));

		}
	
	}
	
	for (int k1=(int)-q/2;k1<=q/2;k1++) {
		exp_k_square_pr[k1+ofs]=exp(-k1*k1/(4*b));
	}


    
    const char* fnames[] = {"nu","alpha","exp_k_square","exp_k1_factor",
                                  "exp_k2_factor","n","m","b","q"};
    mxArray *params=mxCreateStructMatrix(1,1,9,fnames);
    mxSetField(params,0,fnames[0],mxDuplicateArray(nu));
    mxSetField(params,0,fnames[1],alpha);
    mxSetField(params,0,fnames[2],exp_k_square);
    mxSetField(params,0,fnames[3],exp_k1_factor);
    mxSetField(params,0,fnames[4],exp_k2_factor);
    mxSetField(params,0,fnames[5],mxCreateDoubleScalar(n));    
    mxSetField(params,0,fnames[6],mxCreateDoubleScalar(m));
    mxSetField(params,0,fnames[7],mxCreateDoubleScalar(b));
    mxSetField(params,0,fnames[8],mxCreateDoubleScalar(q));
    
    plhs[0] = params;
    
    #ifdef DEBUG
    fclose(fid);
    #endif
    
}
