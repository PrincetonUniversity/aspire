/*

Triangles Scores
[cumulated_scores, scores_histogram] = triangles_scores_mex(N, R_pairs, n_intervals, verbose)

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
For every triangle, the variable
	(m_2nd_best - m_best) / m_2nd_best
which represents the significance of the best configuration's match, turns out to
reveal very interesting information about the triangle and its common lines,
including the rate of indicative common lines among all estimated common lines.
This function wishes to compute the histogram of this variable - the "triangle
score". In addition, for every pair (i,j) it summarizes the scores of all the triangles
{ijk}_k involving (i,j).

Input:
N - number of images
R_pairs - relative rotations matrices (3 X 3 X N-choose-2)
n_intervals - number of intervals in the histogram
verbose - how detailed to print information - unused by now

Output:
cumulated_scores - for every pair ij - the sum of the scores of all the triangles
{ijk} (when k!=i,j) (N-choose-2 X 1)
scores_histogram - histogram of the triangles scores (n_intervals X 1)

Compilation (with inline functions):
mex CFLAGS="\$CFLAGS -std=c99" triangles_scores_mex.c

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
#define CUM_SCORES plhs[0]
#define SCORES_HIST plhs[1]

/* input */
int N = mxGetScalar(prhs[0]);
double *R_pairs = mxGetPr(mxDuplicateArray(prhs[1]));
int n_intervals = mxGetScalar(prhs[2]);
int verbose = mxGetScalar(prhs[3]);

/* initialization */
int N_pairs = N*(N-1)/2;
CUM_SCORES = mxCreateDoubleMatrix(N_pairs,1,mxREAL);
double *cum_scores = mxGetPr(CUM_SCORES);
SCORES_HIST = mxCreateDoubleMatrix(n_intervals,1,mxREAL);
double *scores_hist = mxGetPr(SCORES_HIST);
double h = 1./n_intervals;
int ij,jk,ik;
unsigned int i,j,k,l;
int ALTS[2][4][3];
double s_ij_jk,s_ij_ik,s_ik_jk;
double *Rij, *Rjk, *Rik;
double c[4];
double tmp[9];
int best_i;
double best_val, alt_ij_jk, alt_ik_jk, alt_ij_ik;
double threshold;

/* initialize alternatives */
/* when we find the best J-configuration, we also compare it to the alternative 2nd best one.
 * this comparison is done for every pair in the triplete independently. to make sure that the
 * alternative is indeed different in relation to the pair, we document the differences between
 * the configurations in advance:
 * ALTS(:,best_conf,pair) = the two configurations in which J-sync differs from best_conf in relation to pair */
ALTS[0][0][0]=1; ALTS[0][1][0]=0; ALTS[0][2][0]=0; ALTS[0][3][0]=1; ALTS[1][0][0]=2; ALTS[1][1][0]=3; ALTS[1][2][0]=3; ALTS[1][3][0]=2;
ALTS[0][0][1]=2; ALTS[0][1][1]=2; ALTS[0][2][1]=0; ALTS[0][3][1]=0; ALTS[1][0][1]=3; ALTS[1][1][1]=3; ALTS[1][2][1]=1; ALTS[1][3][1]=1;
ALTS[0][0][2]=1; ALTS[0][1][2]=0; ALTS[0][2][2]=1; ALTS[0][3][2]=0; ALTS[1][0][2]=3; ALTS[1][1][2]=2; ALTS[1][2][2]=3; ALTS[1][3][2]=2;

/* initialize scores */
for (i=0; i<N_pairs; i++) {cum_scores[i] = 0;}
for (i=0; i<n_intervals; i++) {scores_hist[i] = 0;}

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
			
			/* for each triangle side, find the best alternative */
			alt_ij_jk = c[ALTS[0][best_i][0]];
			if (c[ALTS[1][best_i][0]] < alt_ij_jk) {alt_ij_jk = c[ALTS[1][best_i][0]];}
			alt_ik_jk = c[ALTS[0][best_i][1]];
			if (c[ALTS[1][best_i][1]] < alt_ik_jk) {alt_ik_jk = c[ALTS[1][best_i][1]];}
			alt_ij_ik = c[ALTS[0][best_i][2]];
			if (c[ALTS[1][best_i][2]] < alt_ij_ik) {alt_ij_ik = c[ALTS[1][best_i][2]];}
			
			/* compute scores */
			s_ij_jk = 1 - sqrt(best_val/alt_ij_jk);
			s_ik_jk = 1 - sqrt(best_val/alt_ik_jk);
			s_ij_ik = 1 - sqrt(best_val/alt_ij_ik);
			
			/* update cumulated scores */
			cum_scores[ij] += s_ij_jk + s_ij_ik;
			cum_scores[jk] += s_ij_jk + s_ik_jk;
			cum_scores[ik] += s_ik_jk + s_ij_ik;
			
			/* update scores histogram */
			threshold = 0;
			for (l=0; l<n_intervals-1; l++) {
				threshold += h;
				if (s_ij_jk < threshold) {break;}
			}
			scores_hist[l] += 1;
			
			threshold = 0;
			for (l=0; l<n_intervals-1; l++) {
				threshold += h;
				if (s_ik_jk < threshold) {break;}
			}
			scores_hist[l] += 1;
			
			threshold = 0;
			for (l=0; l<n_intervals-1; l++) {
				threshold += h;
				if (s_ij_ik < threshold) {break;}
			}
			scores_hist[l] += 1;
			
		}
	}
}

return;

}
