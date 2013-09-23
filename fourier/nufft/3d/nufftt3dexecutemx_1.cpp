// Function: nufftt3dexecutemx_1
//
// Auxiliry function for nufft_t_3d_execute_1.
// call as nufftt3dexecutevmx(W,precomp);
//
// precomp is the precomputed interpolation tables for the given sampling points.
//
//
// Yoel Shkolnisky, February 2010.

#include <mex.h>
#include <math.h>
#include <string.h>
#include "debug.h"
#include "mexutil.h"

#define for if(0);else for


// Matlab style "mod" function. The result of Matlab's "mod" has the same 
// sign as n. In our case it has to be positive (since n is positive).
inline int matlab_mod(int j, int n) {
    while (j<0) j+=n;
    return j % n;    
}    

void error_msg(const char* msg)
{
    printf("%s: %s",__FILE__,msg);
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
//  SETUPLOG
//	OPENLOG(NULL)
	
    /* if there are no input arguments, print out help information */
    if(nrhs!=2) {
        error_msg("Not enough parameters.\n");
        return;
    }
    const mxArray *W=prhs[0];     // May be complex
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
                                  "exp_k2_factor","exp_k3_factor","n","m","b","q"};

    mxArray *nu= mxGetField(precomp, 0,fnames[0]);
    mxArray *alpha= mxGetField(precomp, 0,fnames[1]);
    mxArray *exp_k_square = mxGetField(precomp, 0,fnames[2]);
    mxArray *exp_k1_factor = mxGetField(precomp, 0,fnames[3]);
    mxArray *exp_k2_factor = mxGetField(precomp, 0,fnames[4]);      
	mxArray *exp_k3_factor = mxGetField(precomp, 0,fnames[5]);      
    GetInteger(mxGetField(precomp, 0,fnames[6]),&n);
    GetInteger(mxGetField(precomp, 0,fnames[7]),&m);
    GetDouble(mxGetField(precomp, 0,fnames[8]),&b);
    GetInteger(mxGetField(precomp, 0,fnames[9]),&q);
	

    const int len_x=mxGetM(alpha);
    const int low_idx_u=(int)ceil((m*n-1.0)/2.0);
    const int high_idx_u=(int)floor((n*m-1.0)/2.0);
   
    mxArray *g=mxCreateDoubleMatrix(len_x, 1, mxCOMPLEX );
    double *g_pr=mxGetPr(g);
    double *g_pi=mxGetPi(g);
    	

    int idx1,idx2,idx3,idx;    
    const int mn=m*n;   
	// Create tables for fast interpolation
        
    double *nu_pr= mxGetPr(nu);
	double *alpha_pr = mxGetPr(alpha);
	double *exp_k_square_pr = mxGetPr(exp_k_square);
	double *exp_k1_factor_pr = mxGetPr(exp_k1_factor);
	double *exp_k2_factor_pr = mxGetPr(exp_k2_factor);
	double *exp_k3_factor_pr = mxGetPr(exp_k3_factor);
    
//	DUMPARR(nu_pr,len_x,3);

	double tmp1,tmp2;
    double Q;
	int    ofs=(int)(q/2);
	double a,c1,c2,c3,d1,d2,d3; // temporary variables to minimize memory access.
	int nu1,nu2,nu3; // temporary variables to minimize memory access.
	double c23,d23;  // temporary variables to save multiplications.
	int idx1base,idx2base,idx3base; // Used to avoid unnecessary additions.
	int ek1,ek2,ek3; // Used to avoid unnecessary additions.
	int eks1,eks2,eks3; // Used to avoid unnecessary additions.
	double gir,gii; // temporary variables to minimize memory access.
	const int q2ofs=-q/2+ofs;

	for (int i=0;i<len_x;i++) {

		a=alpha_pr[i];		
		nu1=(int)nu_pr[i];
		nu2=(int)nu_pr[i+len_x];
		nu3=(int)nu_pr[i+2*len_x];
		idx1base=nu1+low_idx_u;
		idx2base=nu2+low_idx_u;
		idx3base=nu3+low_idx_u;
		ek3=q2ofs*len_x+i;
		eks3=q2ofs;

		gir=0; // Update gir and gii for all k1,k2,k3, and store the results in g_pr and g_pi
		gii=0; // only at the end.

		for (int k3=(int)-q/2;k3<=q/2;k3++) {

			c3=exp_k3_factor_pr[ek3];
			ek3+=len_x;
			d3=exp_k_square_pr[eks3++];

			idx3=matlab_mod(k3+idx3base,mn);
			ek2=q2ofs*len_x+i;
			eks2=q2ofs;

			for (int k2=(int)-q/2;k2<=q/2;k2++) {

				idx2=matlab_mod(k2+idx2base,mn);
				c2=exp_k2_factor_pr[ek2];
				ek2+=len_x;
				d2=exp_k_square_pr[eks2++];

				c23=c2*c3;
				d23=d2*d3;

				ek1=q2ofs*len_x+i;
				eks1=q2ofs;

				for (int k1=(int)-q/2;k1<=q/2;k1++) {

					idx1=matlab_mod(k1+idx1base,mn);					
					c1=exp_k1_factor_pr[ek1];
					ek1+=len_x;
					d1=exp_k_square_pr[eks1++];

					tmp1=d23*d1;
					tmp2=c23*c1;
					Q=a*tmp1*tmp2;

					idx=idx3*mn*mn+idx2*mn+idx1;

					gir=gir+Q*W_pr[idx];
					if (complex_W) {
						gii=gii+Q*W_pi[idx];
					}                
				}
			}
		}
		g_pr[i]=gir;
		if (complex_W) {
			g_pi[i]=gii;
		}                

	}

	plhs[0] = g;    

//	CLOSELOG
}
