// nufftauxmx_v3.cpp - auxiliry function to optimize the performance of 
//                     the polar Fourier transform.
//
// call as nufftauxmx_v3(W,precomp);
//
// precomp is the precomputed interpolation tables for the given sampling points.
//
//
// Yoel Shkolnisky, January 2008.

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
    if(nrhs!=2) {
        error_msg("Not enough parameters.\n");
        return;
    }
    const mxArray *W=prhs[0];     // complex
    const mxArray *precomp=prhs[1];
    
    int n,m;
    double b;
	int q;
    bool complex_W;
        
    if(IsComplexArray(W))
        complex_W=1;
    else if (IsDoubleArray(W))
        complex_W=0;
    else {
        error_msg("W must be a complex or a double array");
        return;
    }
    
    // Get pointers to the data   
    const double *W_pr=mxGetPr(W);   
    double *W_pi;
    
    if (complex_W)
        W_pi=mxGetPi(W);
    else
        W_pi=0;
    
    
    const char* fnames[] = {"nu","alpha","exp_k_square","exp_k1_factor",
                                  "exp_k2_factor","n","m","b","q"};

    mxArray *nu= mxGetField(precomp, 0,fnames[0]);
    mxArray *alpha= mxGetField(precomp, 0,fnames[1]);
    mxArray *exp_k_square = mxGetField(precomp, 0,fnames[2]);
    mxArray *exp_k1_factor = mxGetField(precomp, 0,fnames[3]);
    mxArray *exp_k2_factor = mxGetField(precomp, 0,fnames[4]);      
    GetInteger(mxGetField(precomp, 0,fnames[5]),&n);
    GetInteger(mxGetField(precomp, 0,fnames[6]),&m);
    GetDouble(mxGetField(precomp, 0,fnames[7]),&b);
    GetInteger(mxGetField(precomp, 0,fnames[8]),&q);

    
    const int len_x=mxGetM(alpha);
    const int low_idx_u=(int)ceil((m*n-1.0)/2.0);
    const int high_idx_u=(int)floor((n*m-1.0)/2.0);
   
    mxArray *g=mxCreateDoubleMatrix(len_x, 1, mxCOMPLEX );
    double *g_pr=mxGetPr(g);
    double *g_pi=mxGetPi(g);
       
    int idx1,idx2,idx;    
    const int mn=m*n;   
	// Create tables for fast interpolation
        
    double *nu_pr= mxGetPr(nu);
	double *alpha_pr = mxGetPr(alpha);
	double *exp_k_square_pr = mxGetPr(exp_k_square);
	double *exp_k1_factor_pr = mxGetPr(exp_k1_factor);
	double *exp_k2_factor_pr = mxGetPr(exp_k2_factor);
    
	double tmp1,tmp2;
    double Q;
	int    ofs=(int)(q/2);
     
    for (int k1=(int)-q/2;k1<=q/2;k1++) {
        for (int k2=(int)-q/2;k2<=q/2;k2++) {
			
			tmp1=exp_k_square_pr[k1+ofs]*exp_k_square_pr[k2+ofs];
            
			for (int i=0;i<len_x;i++) {
                
                idx1=(int) nu_pr[i]+k1+low_idx_u;
                while (idx1<0) idx1+=mn;
                idx1=idx1 % mn; // This line together with the previous one implement mod(idx1,mn) with an answer in the range [0,mn-1].
                idx2=(int) nu_pr[i+len_x]+k2+low_idx_u;
                while (idx2<0) idx2+=mn;
                idx2=idx2 % mn;
                idx=idx2*mn+idx1;
                
//                 idx1=((int)nu_pr[i]+k1+low_idx_u+4*mn)%mn;
//                 idx2=((int) nu_pr[i+len_x]+k2+low_idx_u+4*mn)%mn;
//                 idx=idx2*mn+idx1;
                
				tmp2=exp_k1_factor_pr[(k1+ofs)*len_x+i]*exp_k2_factor_pr[(k2+ofs)*len_x+i];
				Q=alpha_pr[i]*tmp1*tmp2;

                g_pr[i]=g_pr[i]+Q*W_pr[idx];
                if (complex_W) {
                    g_pi[i]=g_pi[i]+Q*W_pi[idx];
                }                
            }
        }
    }
    
	plhs[0] = g;    
}
