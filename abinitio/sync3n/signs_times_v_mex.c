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
mex CFLAGS="\$CFLAGS -std=c99" outer_J_sync_mex.c

Written by Ido Greenberg, 2015

*/

#include <math.h>
#include <matrix.h>
#include <time.h>
#include <stdlib.h>
#include <mex.h>


/* from i,j indeces to the common index in the N-choose-2 sized array */
#define PAIR_IDX(N,I,J) ((2*N-I-1)*I/2+J-I-1)

inline void mult_3x3(double *out, double *R1, double *R2) {
/* 3X3 matrices multiplication: out = R1*R2 */
	int i,j;
	for (i=0; i<3; i++) {
		for (j=0;j<3;j++) {
			out[3*j+i] = R1[3*0+i]*R2[3*j+0] + R1[3*1+i]*R2[3*j+1] + R1[3*2+i]*R2[3*j+2];
		}
	}
}

inline void JRJ(double *R) {
/* multiple 3X3 matrix by J from both sizes: R = JRJ */
	R[2]=-R[2];
	R[5]=-R[5];
	R[6]=-R[6];
	R[7]=-R[7];
}

inline double diff_norm_3x3(const double *R1, const double *R2) {
/* difference 2 matrices and return squared norm: ||R1-R2||^2 */
	int i;
	double norm = 0;
	for (i=0; i<9; i++) {norm += (R1[i]-R2[i])*(R1[i]-R2[i]);}
	return norm;
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
/* main */

/* output */
#define V_OUT plhs[0]

/* input */
int N = mxGetScalar(prhs[0]);
double *R_pairs = mxGetPr(mxDuplicateArray(prhs[1]));
double *v = mxGetPr(mxDuplicateArray(prhs[2]));
int N_eigs = mxGetScalar(prhs[3]);
int scores_as_entries = mxGetScalar(prhs[4]);

/* initialization */
int N_pairs = N*(N-1)/2;
V_OUT = mxCreateDoubleMatrix(N_eigs*N_pairs,1,mxREAL);
double *v2 = mxGetPr(V_OUT);
int ij,jk,ik;
unsigned int i,j,k;
double s_ij_jk,s_ij_ik,s_ik_jk;
int signs_confs[4][3], ALTS[2][4][3];
double *Rij, *Rjk, *Rik;
double c[4];
double tmp[9];
int best_i;
double best_val, alt_ij_jk, alt_ik_jk, alt_ij_ik;
int eig;

/* initialize signs configurations and alternatives */
for(i=0; i<4; i++) { for(j=0; j<3; j++) { signs_confs[i][j]=1; } }
signs_confs[2-1][1-1]=-1; signs_confs[2-1][3-1]=-1;
signs_confs[3-1][1-1]=-1; signs_confs[3-1][2-1]=-1;
signs_confs[4-1][2-1]=-1; signs_confs[4-1][3-1]=-1;
/* initialize alternatives */
/* when we find the best J-configuration, we also compare it to the alternative 2nd best one.
 * this comparison is done for every pair in the triplete independently. to make sure that the
 * alternative is indeed different in relation to the pair, we document the differences between
 * the configurations in advance:
 * ALTS(:,best_conf,pair) = the two configurations in which J-sync differs from best_conf in relation to pair */
ALTS[0][0][0]=1; ALTS[0][1][0]=0; ALTS[0][2][0]=0; ALTS[0][3][0]=1; ALTS[1][0][0]=2; ALTS[1][1][0]=3; ALTS[1][2][0]=3; ALTS[1][3][0]=2;
ALTS[0][0][1]=2; ALTS[0][1][1]=2; ALTS[0][2][1]=0; ALTS[0][3][1]=0; ALTS[1][0][1]=3; ALTS[1][1][1]=3; ALTS[1][2][1]=1; ALTS[1][3][1]=1;
ALTS[0][0][2]=1; ALTS[0][1][2]=0; ALTS[0][2][2]=1; ALTS[0][3][2]=0; ALTS[1][0][2]=3; ALTS[1][1][2]=2; ALTS[1][2][2]=3; ALTS[1][3][2]=2;

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

            /* compute configurations matches scores */
            Rij = R_pairs + 9*ij;
            Rjk = R_pairs + 9*jk;
            Rik = R_pairs + 9*ik;

            mult_3x3(tmp,Rij,Rjk);
            c[0] = diff_norm_3x3(tmp,Rik);

            JRJ(Rij);
            mult_3x3(tmp,Rij,Rjk);
            c[1] = diff_norm_3x3(tmp,Rik);
            JRJ(Rij);

            JRJ(Rjk);
            mult_3x3(tmp,Rij,Rjk);
            c[2] = diff_norm_3x3(tmp,Rik);
            JRJ(Rjk);

            JRJ(Rik);
            mult_3x3(tmp,Rij,Rjk);
            c[3] = diff_norm_3x3(tmp,Rik);
            JRJ(Rik);

            /* find best match */
            best_i=0; best_val=c[0];
            if (c[1]<best_val) {best_i=1; best_val=c[1];}
            if (c[2]<best_val) {best_i=2; best_val=c[2];}
            if (c[3]<best_val) {best_i=3; best_val=c[3];}

            if (scores_as_entries != 0) {
                /* for each triangle side, find the best alternative */
                alt_ij_jk = c[ALTS[0][best_i][0]];
                if (c[ALTS[1][best_i][0]] < alt_ij_jk) {alt_ij_jk = c[ALTS[1][best_i][0]];}
                alt_ik_jk = c[ALTS[0][best_i][1]];
                if (c[ALTS[1][best_i][1]] < alt_ik_jk) {alt_ik_jk = c[ALTS[1][best_i][1]];}
                alt_ij_ik = c[ALTS[0][best_i][2]];
                if (c[ALTS[1][best_i][2]] < alt_ij_ik) {alt_ij_ik = c[ALTS[1][best_i][2]];}

                /* compute triangles entries */
                s_ij_jk = signs_confs[best_i][0] * (1 - sqrt(best_val/alt_ij_jk));
                s_ik_jk = signs_confs[best_i][1] * (1 - sqrt(best_val/alt_ik_jk));
                s_ij_ik = signs_confs[best_i][2] * (1 - sqrt(best_val/alt_ij_ik));
            }
            else {
                /* set triangles entries */
                s_ij_jk = signs_confs[best_i][0];
                s_ik_jk = signs_confs[best_i][1];
                s_ij_ik = signs_confs[best_i][2];
            }


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
