// nufftauxmx_v2.cpp - auxiliry function to optimize the performance of pft.m
//
// call as nufftauxmx_v2(x,W,nu,n,m,b,q);
//
// This function is an optimized version of nufftauxmx. It pre-computes the interpolation factors to 
// avoid unnecessary exponent evaluations. These are stored in the tables
//	alpha, exp_k_square, exp_k1_factor, exp_k2_factor.
//
//
// Yoel Shkolnisky, May 2007.



//#define DEBUG

#include <mex.h>
#include <math.h>
#include <string.h>
#include <stdio.h>

#define for if(0);else for

/*
 * MATLAB utlity routines
 */


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

void error_msg(const char* msg)
{
    printf("%s: %s",__FILE__,msg);
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
    /* if there are no input arguments, print out help information */
    if(nrhs!=7) {
        error_msg("Not enough parameters.\n");
        return;
    }
    
    #ifdef DEBUG
    FILE *fid=fopen("debug.txt","w");
    #endif
    
    const double PI=atan(1.0)*4;
    
    const mxArray *x=prhs[0];     // double
    const mxArray *W=prhs[1];     // complex
    const mxArray *nu=prhs[2];    // double
    
    int n,m;
    double b;
	int q;
    bool complex_W;
    
    GetInteger(prhs[3],&n);
    GetInteger(prhs[4],&m);
    GetDouble(prhs[5],&b);
    GetInteger(prhs[6],&q);
    
    if(!IsDoubleArray(x)) {
        error_msg("x must be a double array");
        return;
    }
    
    if(IsComplexArray(W))
        complex_W=1;
    else if (IsDoubleArray(W))
        complex_W=0;
    else {
        error_msg("W must be a complex or a double array");
        return;
    }
    
    if(!IsDoubleArray(nu)) {
        error_msg("nu must be a double array");
        return;
    }
    
    // Get pointers to the data
    
    const double *x_pr=mxGetPr(x);
    const double *W_pr=mxGetPr(W);
    const double *nu_pr=mxGetPr(nu);
    
    double *W_pi;
    
    if (complex_W)
        W_pi=mxGetPi(W);
    else
        W_pi=0;
    
    
    const int len_x=mxGetM(x);
    const int low_idx_u=(int)ceil((m*n-1.0)/2.0);
    const int high_idx_u=(int)floor((n*m-1.0)/2.0);
    
    #ifdef DEBUG
    fprintf(fid,"low_idx_u=%d\n",low_idx_u);
    fprintf(fid,"m=%d\n",m);
    fprintf(fid,"n=%d\n",n);
    fprintf(fid,"len_x=%d\n",len_x);
    fprintf(fid,"expr=%7.4f\n",ceil((m*n-1.0)/2.0));
    #endif
    
    mxArray *g=mxCreateDoubleMatrix(len_x, 1, mxCOMPLEX );
    double *g_pr=mxGetPr(g);
    double *g_pi=mxGetPi(g);
       
    int idx1,idx2,idx;    
    double t1,t2,Q;
    const int mn=m*n;
    const double c1=4*b*PI;
    const double c2=mn/(2*PI);
    
/*
    #ifdef DEBUG
    fprintf(fid,"nu= %d x %d\n",mxGetM(nu),mxGetN(nu));
    for (int i=0; i<len_x; i++) {
        fprintf(fid,"%d   %d\n",nu_pr[i],nu_pr[i+len_x]);
    }
    #endif
 */
	// Create tables for fast interpolation
	double *alpha;
	double *exp_k_square;
	double *exp_k1_factor;
	double *exp_k2_factor;
	double tmp,tmp1,tmp2,tmp3;
	int    ofs=(int)(q/2);

    if ((alpha=(double*) malloc(len_x*sizeof(double)))==NULL) {
        error_msg("Failed to allocate alpha...aborting\n");
        return;
    }
    
	if ((exp_k_square=(double*) malloc((q+1)*sizeof(double)))==NULL) {
        error_msg("Failed to allocate exp_k_square...aborting\n");
        return;
    }
    
	if ((exp_k1_factor=(double*) malloc((q+1)*len_x*sizeof(double)))==NULL) {
        error_msg("Failed to allocate exp_k1_factor...aborting\n");        
        return;
    }
    
	if ((exp_k2_factor=(double*) malloc((q+1)*len_x*sizeof(double)))==NULL) {
        error_msg("Failed to allocate exp_k2_factor...aborting\n");
        return;
    }
	
    
    for (int i=0;i<len_x;i++) {

		tmp1=(x_pr[i]*c2)*(x_pr[i]*c2)+(x_pr[i+len_x]*c2)*(x_pr[i+len_x]*c2);
		tmp2=nu_pr[i]*nu_pr[i]+nu_pr[i+len_x]*nu_pr[i+len_x];
		tmp3=-2*x_pr[i]*c2*nu_pr[i]-2*x_pr[i+len_x]*c2*nu_pr[i+len_x];
		tmp=tmp1+tmp2+tmp3;

		alpha[i]=exp(-tmp/(4*b))/c1;

		for (int k=(int)-q/2;k<=q/2;k++) {
			exp_k1_factor[(k+ofs)*len_x+i]=exp(-(-k*(2*x_pr[i]*c2-2*nu_pr[i]))/(4*b));
			exp_k2_factor[(k+ofs)*len_x+i]=exp(-(-k*(2*x_pr[i+len_x]*c2-2*nu_pr[i+len_x]))/(4*b));

		}
	
	}
	
	for (int k1=(int)-q/2;k1<=q/2;k1++) {
		exp_k_square[k1+ofs]=exp(-k1*k1/(4*b));
	}

	
    for (int k1=(int)-q/2;k1<=q/2;k1++) {
        for (int k2=(int)-q/2;k2<=q/2;k2++) {
			
			tmp1=exp_k_square[k1+ofs]*exp_k_square[k2+ofs];
            
			for (int i=0;i<len_x;i++) {
 /*               
                idx1=(int) nu_pr[i]+k1+low_idx_u;
                idx1=idx1+((int)(abs(idx1/mn)+1))*mn;
                idx1=idx1 % mn; // This line together with the previous one implement mod(idx1,mn) with an answer in the range [0,mn-1].
                idx2=(int) nu_pr[i+len_x]+k2+low_idx_u;
                idx2=idx2+((int)(abs(idx2/mn)+1))*mn;
                idx2=idx2 % mn;
                idx=idx2*mn+idx1;
*/
                idx1=(int) nu_pr[i]+k1+low_idx_u;
				while (idx1<0) idx1+=mn;
                idx1=idx1 % mn; // This line together with the previous one implement mod(idx1,mn) with an answer in the range [0,mn-1].
                idx2=(int) nu_pr[i+len_x]+k2+low_idx_u;
				while (idx2<0) idx2+=mn;
                idx2=idx2 % mn;
                idx=idx2*mn+idx1;


				tmp2=exp_k1_factor[(k1+ofs)*len_x+i]*exp_k2_factor[(k2+ofs)*len_x+i];
				Q=alpha[i]*tmp1*tmp2;

                g_pr[i]=g_pr[i]+Q*W_pr[idx];
                if (complex_W) {
                    g_pi[i]=g_pi[i]+Q*W_pi[idx];
                }
                
                #ifdef DEBUG
                fprintf(fid,"k1=%d k2=%d i=%d j=%d  tmp=%7.4e  Q=%7.4e   g=%7.4e\n",k1,k2,i,idx,tmp,Q,g_pr[i]);
                #endif
            }
        }
    }
    
	free(alpha);
	free(exp_k_square);
	free(exp_k1_factor);
	free(exp_k2_factor);
    
	plhs[0] = g;

    #ifdef DEBUG
    fclose(fid);
    #endif
    
}
