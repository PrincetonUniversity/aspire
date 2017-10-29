// bfsnn (Brute force nearest neighbor search)
//
//
//
// Each column correponds to a point. Each row is a coordinate.
//
// Yoel Shkolnisky, May 2007.

//#define DEBUG

#include <mex.h>
#include <math.h>
#include <string.h>


#define IMAGE_NAME "ksden"

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


void usage()
{
	printf("Density estimation using a Gaussian kernel.\n");
	printf("d=ksden(X,sigma)\n");
	printf("\n");
	printf("   X      M by N dataset. N data points. M coordinates to each data point\n");
	printf("   sigma  Width of the Guassian kernel\n");
	printf("\n");
	printf("Returned values:\n");	
	printf("   d      Vector of length N with entry i being the density estimation at point i.\n");
}



void error_msg(const char* msg)
{
	printf("%s: %s",IMAGE_NAME,msg);
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{   
	
    /* if there are no input arguments, print out help information */
    if(nrhs==0) {
        usage();
        return;
    }       
    
    if(nrhs < 2) {
        error_msg("one or more input parameters missing.\n");
        return;
    }
	
    if(nrhs > 2) {
        error_msg("too many input parameters.\n");
        return;
    }
	
    const mxArray *dataset = prhs[0];    
    if(!mxIsDouble(dataset) | mxIsSparse(dataset) | mxIsComplex(dataset)) {
        error_msg("dataset must be a real full matrix of doubles.\n");
        return;
    }
	
    int M=mxGetM(dataset); // Each data point has M coordinates
	int N=mxGetN(dataset); // N data points
	double *dataset_pr=mxGetPr(dataset);
	
	double sigma;
	GetDouble(prhs[1],&sigma);
	
	mxArray *density    = mxCreateDoubleMatrix(N , 1, mxREAL);
	double  *density_pr = mxGetPr(density);
	
	double *row_dists=(double*)malloc(N*sizeof(double));
	
	double sum;
	double temp;
	double dij;


	for (int i=0; i<N; i++) {				
		// compute all distances for point i		
		for (int j=0; j<i; j++) {

			sum=0.0;			

			for(int k=0; k < M; k++) {
				temp = dataset_pr[i*M+k]-dataset_pr[j*M+k];
				sum += temp*temp;
			}
			
			dij=exp(-sum/sigma);
			density_pr[i]+=dij;			
			density_pr[j]+=dij;
		}					
	}
	
	
	// Sort each column of distances and idx.
	
	free(row_dists);
	
	plhs[0] = density;
}
