// nufftauxmx_t3.cpp - auxiliry function to optimize the performance of 
//                     nufft_t3_execute
//
// call as nufftauxmx_t3(x,precomp);
//
// Yoel Shkolnisky, January 2009.

#define IMAGE_NAME "nufftauxmx_t3"

#include <mex.h>
#include <math.h>
#include <string.h>
#include "../common/mexutil.h"

#define for if(0);else for

// Increment modulo n.
#define MODINC(j,n)\
     (j)++;\
      if ((j)==(n))\
           (j)=0;\


// Matlab style "mod" function. The result of Matlab's "mod" has the same 
// sign as n. In our case it has to be positive (since n is positive).
inline int matlab_mod(int j, int n) {
    while (j<0) j+=n;
    return j % n;    
}    

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{    
    if(nrhs!=2) {
        mexErrMsgTxt("Not enough parameters.\n");
     }
    const mxArray *x=prhs[0];
    const mxArray *precomp=prhs[1];
    
    // Check the dimensions of x.
    if (mxGetNumberOfDimensions(x)!=2) {
        mexErrMsgTxt("x must be a mx3 array");
    }
    if (mxGetN(x)!=3) {
        mexErrMsgTxt("x must be a mx3 array");
    }
        
    // x contains points in 3D so there is no imaginary part.
    const double *x_pr=mxGetPr(x);   
    const int len_x=mxGetM(x);
    
    // Read the fields from the precomp structure.
    const char* fnames[] = {"W","b","m","q","n"}; // Field names.
    double b;
    int m,n,q;
    
    mxArray *W= mxGetField(precomp, 0,fnames[0]);
    GetDouble(mxGetField(precomp, 0,fnames[1]),&b);
    GetInteger(mxGetField(precomp, 0,fnames[2]),&m);
    GetInteger(mxGetField(precomp, 0,fnames[3]),&q);
    GetInteger(mxGetField(precomp, 0,fnames[4]),&n);      
         
    double *W_pr=mxGetPr(W);
    double *W_pi=mxGetPi(W);
       
    // Allocate output array.
    mxArray *g=mxCreateDoubleMatrix(len_x, 1, mxCOMPLEX );
    double  *g_pr=mxGetPr(g);
    double  *g_pi=mxGetPi(g);
        
    // Precompute recurring constants 
    const int mn=m*n;
    const int low_idx_u=ceil((mn-1)/2.0);
    const double PI=atan(1.0)*4;
    const double c1=pow(2*sqrt(b*PI),3);
    const int c2=-q/2+low_idx_u;
    
    
    for (int k=0; k<len_x ; k++) {
        // Take current row from x.
        double x0=x_pr[k];
        double x1=x_pr[k+len_x];
        double x2=x_pr[k+2*len_x];
                
        int nu0=round(x0*mn/(2*PI));
        int nu1=round(x1*mn/(2*PI));
        int nu2=round(x2*mn/(2*PI));
                
        double xx0=x0*mn/(2*PI)-nu0;
        double xx1=x1*mn/(2*PI)-nu1;
        double xx2=x2*mn/(2*PI)-nu2;
           
        int j3=matlab_mod(nu2+c2,mn);
        
        for (int k3=-q/2; k3<=q/2; k3++) {               
            int j2=matlab_mod(nu1+c2,mn);
            
            for (int k2=-q/2; k2<=q/2; k2++) {               
                int jj=(j3*mn+j2)*mn; 
                int j1=matlab_mod(nu0+c2,mn);
                int j = jj+j1;
                
                for (int k1=-q/2; k1<=q/2; k1++) {                    
                    double tmp1=(xx0-k1);
                    double tmp2=(xx1-k2);
                    double tmp3=(xx2-k3);                    
                    double tmp= -(tmp1*tmp1+tmp2*tmp2+tmp3*tmp3)/(4*b);
                    double Q=exp(tmp);
                                        
                    g_pr[k]=g_pr[k]+Q*W_pr[j]/c1;
                    g_pi[k]=g_pi[k]+Q*W_pi[j]/c1;
                                         
                    j1++; j++;
                    
                    if(j1==mn){
                        j1=0;
                        j=jj;
                    }
                }             
                MODINC(j2,mn);
            }                                    
            MODINC(j3,mn);
        }                
    }
    
    plhs[0] = g;
}


/* Reference Matlab code:
 *
 * % g=zeros(len_x,1);
 * % mn=m*n;
 * % c1=(2*sqrt(b*pi))^3;
 * % low_idx_u=ceil((m*n-1)/2);
 * % 
 * % for k=1:len_x
 * %     
 * %     nu=round(x(k,:)*m*n/(2*pi));
 * %     xx=x(k,:).*m.*n/(2*pi)-nu;
 * %     
 * %     j3=mod(nu(3)-q/2+low_idx_u,mn)+1;
 * %     for k3=-q/2:q/2
 * %         j2=mod(nu(2)-q/2+low_idx_u,mn)+1;
 * %         for k2=-q/2:q/2
 * %             jj=(j3-1)*(mn)^2 + (j2-1)*mn;
 * %             j1=mod(nu(1)-q/2+low_idx_u,mn)+1;
 * %             j=jj+j1;
 * %             
 * %             for k1=-q/2:q/2
 * %                 tmp=-((xx(1)-k1)^2+(xx(2)-k2)^2+(xx(3)-k3)^2)/(4*b);
 * %                 Q=exp(tmp);
 * %                 g(k)=g(k)+Q.*W(j)/c1;
 * %                 
 * %                 j1=j1+1;
 * %                 j=j+1;
 * %                 if j1==mn+1
 * %                     j1=1;
 * %                     j=jj+1;
 * %                 end
 * %             end
 * %             
 * %             j2=j2+1;
 * %             if j2==mn+1
 * %                 j2=1;
 * %             end
 * %         end
 * %                 
 * %         j3=j3+1;
 * %         if j3==m*n+1
 * %             j3=1;
 * %         end
 * %     end
 * % end
 */


