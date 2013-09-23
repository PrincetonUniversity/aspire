/*
Gridding loop for the 3D NUFFT (reference implementation).

This is a rather slow (but readable) implementation obtained by a direct translation of the Matlab code attached below. 

Example:
	tau=nufft3dauxmx_ref(n,M,m,q,mu,Px,Py,Pz,alpha);

Yoel Shkolnisky, January 2010.
*/

#include <mex.h>
#include <math.h>
#include <string.h>
#include "mexutil.h"

#define for if(0);else for

inline int matlab_mod(int j, int n) {
	while (j<0) j+=n;
	return j % n;    
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	// Parameters: n,M,m,q,mu,Px,Py,Pz,alpha
	// Returns: tau

	int n,M,m,q;
	double *mu,*Px,*Py,*Pz,*alpha_pr,*alpha_pi;	
	double *tau_pr, *tau_pi; // output variables

	mxComplexity complexFlag; // Is the input array alpha complex?
	int dims[3];
	int j1,j2,j3,idx1,idx2,idx3,idx,k,offset;
	double w1,w2,w12;

	if (nrhs!=9) {
		mexErrMsgTxt("9 parameters are required: n,M,m,q,mu,Px,Py,Pz,alpha");
	}

	// Get Input parameters
	GetInteger(prhs[0],&n);
	GetInteger(prhs[1],&M);
	GetInteger(prhs[2],&m);
	GetInteger(prhs[3],&q);
	mu=mxGetPr(prhs[4]);
	Px=mxGetPr(prhs[5]);
	Py=mxGetPr(prhs[6]);
	Pz=mxGetPr(prhs[7]);
	alpha_pr=mxGetPr(prhs[8]);

	complexFlag=mxREAL;
	alpha_pi=NULL;
	if (mxIsComplex(prhs[8])) {
		alpha_pi=mxGetPi(prhs[8]);
		complexFlag=mxCOMPLEX;
	}

	// Allocate output variable
	dims[0]=m*M; dims[1]=m*M; dims[2]=m*M;
	plhs[0]=mxCreateNumericArray(3,dims,mxDOUBLE_CLASS,complexFlag);
	tau_pr=mxGetPr(plhs[0]);
	tau_pi=mxGetPi(plhs[0]);


	// Gridding loop
	offset=(int)ceil((m*M-1.0)/2.0);
	for (k=0; k<n; k++) {
		for (j1=-q/2; j1<=q/2; j1++) {
			idx1=(int)mu[k]+j1;
			idx1=matlab_mod(idx1+offset,m*M);
			w1=Px[n*(j1+q/2)+k];
			for (j2=-q/2; j2<=q/2; j2++) {
				idx2=(int)mu[n+k]+j2;
				idx2=matlab_mod(idx2+offset,m*M);
				w2=Py[n*(j2+q/2)+k];
				w12=w1*w2;

				for (j3=-q/2; j3<=q/2; j3++) {
					idx3=(int)mu[2*n+k]+j3;
					idx3=matlab_mod(idx3+offset,m*M);
					idx=idx3*(m*M)*(m*M)+idx2*m*M+idx1;
					tau_pr[idx]+=w12*Pz[n*(j3+q/2)+k]*alpha_pr[k];
					if (complexFlag) {
						tau_pi[idx]+=w12*Pz[n*(j3+q/2)+k]*alpha_pi[k];
					}
				}
			}
		}    

	}
}

// Matlab reference code for the above main loop:
//for j1=-q/2:q/2
//    idx1=mu(:,1)+j1;
//    idx1=mod(idx1+offset1,m*M);
//    w1=Px(:,j1+q/2+1);
//    for j2=-q/2:q/2
//        idx2=mu(:,2)+j2;
//        idx2=mod(idx2+offset1,m*M);
//        w2=Py(:,j2+q/2+1);
//        w12=w1.*w2;
//        for j3=-q/2:q/2
//            idx3=mu(:,3)+j3;
//            idx3=mod(idx3+offset1,m*M);
//            idx=idx3.*(m*M)^2+idx2.*m*M+idx1+1; % convert to linear index
//            tau(idx(:))=tau(idx(:))+...
//                w12.*Pz(:,j3+q/2+1).*alpha(:);
//        end
//    end
//end
