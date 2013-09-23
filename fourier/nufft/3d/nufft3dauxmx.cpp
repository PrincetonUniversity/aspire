/*
MEX version of the gridding loop for the 3D NUFFT

Optimized version of nufftdauxmx_ref.cpp.
The iplemented optimizations are:
1) Move the varaible k to the outermost loop.
2) Remove unnecessary access to memory. For example, access mu only from the outer loop.
3) The arrays Px Py and Pz are transposed compared to nufftdauxmx_ref, to access memory in columns.
4) The loops for j1,j2,j3 now go from 0 to q instead of -q/2 to q/2.
5) Instad of computing the indices to Px, Py, Pz each time, we maintain the indices idxPx,idxPy,idxPz, and just update them each iteration.

Example:
	tau=nufft3dauxmx(n,M,m,q,mu,Px.',Py.',Pz.',alpha);
Note the transpose of Px, Py, and Pz compared to nufft3dauxmx_ref.

Yoel Shkolnisky, January 2010.
*/

#include <mex.h>
#include <math.h>
#include <string.h>
#include "../common/mexutil.h"

#define for if(0);else for


inline int matlab_mod(int j, int n) {
	while (j<0) j+=n;
	return j % n;    
//	The above code is no slower than "return (j+4*n) % n".
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	// Parameters: n,M,m,q,mu,Px,Py,Pz,alpha
	// Returns: tau

	int n,M,m,q;
	double *mu,*Px,*Py,*Pz,*alpha_pr,*alpha_pi;	
	double *tau_pr, *tau_pi; // output variable

	mxComplexity complexFlag; // Is the input array alpha complex?
	int dims[3];
	
	int k,j1,j2,j3;
	// We regrid each point omega(k) to (q+1)^3 regular points around it.
	// The regular grid points are the q+1 points closest to round(m*omega(k)).
	// The variables startinggridpointx, startinggridpointy, startinggridpointz hold 
	// the coordinates of first grid point of the regular grid around omega(k), in the x,y, 
	// and z directions, respectively. The coordinates are given in terms of the padded volume 
	// tau of size (m*M)^3.
	int startinggridpointx,startinggridpointy,startinggridpointz; 
	// Current grid point processed. "gp" stands for "grid point". (gpx,gpy,gpz) is the current 
	// grid point of the (q+1)^3 regular cube around omega(k), at which we compute the conribution 
	// of the point omega(k).
	int gpx,gpy,gpz;
	// Weights to be used at the point (gpx,gpy,gpz). idxPx is the weight used to compute the contibution
	// of the point omega(k) at the point gpx. idxPy and idxPz have similar roles.
	int idxPx,idxPy,idxPz;	
	int gpxm,gpym,gpzm; // The points gpx,gpy,gpz modulu m*M.
	int idx; // The linear index of the point (gpxm,gpym,gpzm) in an array of size (m*M)^3.
	// Absolute value of the lowest index in a zero-centered array. For example, "2" if we have a 5 elements 
	// array, since the indices are -2,1,0,1,2.  This variable is used to translate indices of zero-centered 
	// arrays to indices of zero-based arrays.
	int offset; 
	double w,w1,w2,w12; // Temporary variables for computing the gridding weights.
	int mM,mM2; // To avoid repeated computations of m*M and (m*M)^2.
	double akr,aki; // Aviod accessing alpha_pr and alpha_pi for each j1,j2,j3. Access it only once.


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

	// Gridding loop.
	mM=m*M;
	mM2=mM*mM;
	offset=(int)ceil((m*M-1.0)/2.0);
	for (k=0; k<n; k++) {
		startinggridpointx=(int)mu[k]-q/2;
		startinggridpointy=(int)mu[n+k]-q/2;
		startinggridpointz=(int)mu[2*n+k]-q/2;
		gpx=startinggridpointx+offset;
		idxPx=k*(q+1);

		akr=alpha_pr[k];
		if (complexFlag) {
			aki=alpha_pi[k];
		}

		for (j1=0; j1<=q; j1++) {
			gpxm=matlab_mod(gpx,mM);			
			w1=Px[idxPx];
			gpx++;
			idxPx++;
			gpy=startinggridpointy+offset;			
			idxPy=k*(q+1);
			
			for (j2=0; j2<=q; j2++) {
				gpym=matlab_mod(gpy,mM);				
				w2=Py[idxPy];
				gpy++;
				idxPy++;
				w12=w1*w2;
				gpz=startinggridpointz+offset;	
				idxPz=k*(q+1);

				for (j3=0; j3<=q; j3++) {
					gpzm=matlab_mod(gpz,mM);					
					idx=gpzm*mM2+gpym*mM+gpxm; // There is no point in doing the computation incrementaly.
					w=w12*Pz[idxPz];
					tau_pr[idx]+=w*akr;
					if (complexFlag) {
						tau_pi[idx]+=w*aki;
					}
					gpz++;
					idxPz++;
				}
			}
		}    
	}

}
