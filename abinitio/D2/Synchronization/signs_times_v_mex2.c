/*

Outer J Sync
[v, eigen_value, n_iterations] = outer_J_sync_mex(N, R_pairs, v0, verbose)

Background:
Common lines between images derive the relative rotations between the images (Rij)
only up to a reflection, i.e. for each pair ij we have either Rij OR J*Rij*J.
This code synchronizes the reflections through the relation Rij*Rjk*Rki=I,
which should hold for every triangle ijk. In fact, we just try variant configurations
of J-multiplications of the rotations Rij/Rjk/Rki, and choose the one that minimizes
m=match=||Rij*Rjk*Rki-I||.
This way, for every such triangle we find the most likely signs configuration, in which
1 represents the same J presence in both rotations, and -1 means that only one of the
rotations should be multiplied by J.
Since we have a synchronization for every pair of pairs ij,jk, we can represent the
data by N-choose-2 on N-choose-2 matrix. The global synchronization of the rotations
is given by its first N-choose-2 sized eigenvector, in which 1 means Rij is correct,
and -1 means that it should be multiplied by J (or vice versa - the synchronization is
inherently only up to a global reflection).

This function's goal:
To avoid a huge memory usage (~N^3), we do not keep the entries of the signs matrix.
Instead, we use the power method for the eigenvector and re-compute all the entries
every iteration to apply the multiplication.
This function implements the computation of the entries of the signs matrix, and the
multiplication of this matrix by the eigenvector-candidate.
It is called by every iteration of the power method.

Input:
N - number of projections
R_pairs - relative rotations matrices (3 X 3 X N-choose-2)
v - power method's current eigenvectors-candidates (N-choose-2 X N_eigs)
N_eigs - number of eigenvectors computed
scores_as_entries - 0 for +-1 matrix entries, non 0 for weighted entries

Output:
V_OUT - signs_matrix*v multiplication (N-choose-2*N_eigs X 1,
 which will be reshaped as N-choose-2 X N_eigs)

Compilation (with inline functions):
mex CFLAGS="\$CFLAGS -std=c99" signs_times_v_mex.c

Written by Ido Greenberg, 2015
Adapted for D2 by Eitan Rosen, 2016

*/

#include <math.h>
#include <matrix.h>
#include <time.h>
#include <stdlib.h>
#include <mex.h>
#include <tmwtypes.h>


/* from i,j indeces to the common index in the N-choose-2 sized array */
#define PAIR_IDX(N,I,J) ((2*N-I-1)*I/2+J-I-1)
#define TRIP_IDX(N,I,J,K) (N*(N-1)*(N-2)/6-(N-I-3)*(N-I-2)*(N-I-1)/6-(N-J-1)*(N-J)/2+K-J-1);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
/* main */

/* output */
#define V_OUT plhs[0]  

/* input */    
unsigned long long N = mxGetScalar(prhs[0]);
double *best_vec = mxGetPr(mxDuplicateArray(prhs[1]));
double *v = mxGetPr(mxDuplicateArray(prhs[2]));
int N_eigs = mxGetScalar(prhs[3]);

/* initialization */
unsigned long long N_pairs = N*(N-1)/2;
unsigned long long ij,jk,ik,ijk;
unsigned long long i,j,k;
int s_ij_jk,s_ij_ik,s_ik_jk;
int signs_confs[4][3];
unsigned int best_i;
int eig;
double *v2;

V_OUT = mxCreateDoubleMatrix(N_eigs*N_pairs,1,mxREAL);
v2 = mxGetPr(V_OUT);

/* initialize signs configurations and alternatives */
for(i=0; i<4; i++) { for(j=0; j<3; j++) { signs_confs[i][j]=1; } }
signs_confs[2-1][1-1]=-1; signs_confs[2-1][3-1]=-1;
signs_confs[3-1][1-1]=-1; signs_confs[3-1][2-1]=-1;
signs_confs[4-1][2-1]=-1; signs_confs[4-1][3-1]=-1;

/* initialize v2 */
for (i=0; i<N_eigs*N_pairs; i++) {v2[i]=0;}

/* loop all triangles */
for (i=0; i<(N-2); i++) {
    for (j=i+1; j<(N-1); j++) {
        ij = PAIR_IDX(N,i,j);
        jk = PAIR_IDX(N,j,j);
        ik = ij;
        for (k=j+1; k<N; k++) {
            jk++;
            ik++;
            
            ijk=TRIP_IDX(N,i,j,k);
            best_i=(unsigned int) best_vec[ijk];                
            s_ij_jk = signs_confs[best_i][0];
            s_ik_jk = signs_confs[best_i][1];
            s_ij_ik = signs_confs[best_i][2];
            
            /* update multiplication */
            for (eig=0; eig<N_eigs; eig++) {
                v2[ij+N_pairs*eig] += s_ij_jk*v[jk+N_pairs*eig] + s_ij_ik*v[ik+N_pairs*eig];
                v2[jk+N_pairs*eig] += s_ij_jk*v[ij+N_pairs*eig] + s_ik_jk*v[ik+N_pairs*eig];
                v2[ik+N_pairs*eig] += s_ij_ik*v[ij+N_pairs*eig] + s_ik_jk*v[jk+N_pairs*eig];
            }
        }
    }
}

return;

}
