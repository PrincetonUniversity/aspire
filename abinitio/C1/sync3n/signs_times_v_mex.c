/*

Signs x V Multiplication
[v, cores_used] = signs_times_v_mex...
    (N, R_pairs, v0, n_eigs, scores_as_entries, pairs_scores, n_threads);

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
scores_as_entries - 0 for +-1 matrix entries, 1 for triangles scores entries,
2 for pair-wise entries
pairs_scores - if scores_as_entries==2, then these scores used to weight the
entries of the triplets (N-choose-2 array)
n_threads - how many cores to use (if available)

Output:
V_OUT - signs_matrix*v multiplication (N-choose-2*N_eigs X 1,
 which will be reshaped as N-choose-2 X N_eigs)
cores_used - how many cores were used in practice

Compilation (with inline functions):
mex CFLAGS="\$CFLAGS -std=c99" signs_times_v_mex.c

Written by Ido Greenberg, 2016

*/

#include <math.h>
#include <matrix.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <mex.h>
#include <pthread.h>


/* from i,j indeces to the common index in the N-choose-2 sized array */
#define PAIR_IDX(N,I,J) ((2*N-I-1)*I/2+J-I-1)

struct ThreadData {
	unsigned int start, stop;
	unsigned long N;
	double *R_pairs;
	double *v;
	unsigned int N_eigs;
	int scores_as_entries;
	double * pairs_scores;
	double *v2;
	int (*ALTS)[4][3];
	int (*signs_confs)[3];
};

inline void mult_3x3(double *out, double *R1, double *R2) {
/* 3X3 matrices multiplication: out = R1*R2 */
	int i,j;
	for (i=0; i<3; i++) {
		for (j=0;j<3;j++) {
			out[3*j+i] = R1[3*0+i]*R2[3*j+0] + R1[3*1+i]*R2[3*j+1] + R1[3*2+i]*R2[3*j+2];
		}
	}
}

inline void JRJ(double *R, double *A) {
/* multiple 3X3 matrix by J from both sizes: A = JRJ */
	A[0]=R[0];
	A[1]=R[1];
	A[2]=-R[2];
	A[3]=R[3];
	A[4]=R[4];
	A[5]=-R[5];
	A[6]=-R[6];
	A[7]=-R[7];
	A[8]=R[8];
}

inline double diff_norm_3x3(const double *R1, const double *R2) {
/* difference 2 matrices and return squared norm: ||R1-R2||^2 */
	int i;
	double norm = 0;
	for (i=0; i<9; i++) {norm += (R1[i]-R2[i])*(R1[i]-R2[i]);}
	return norm;
}

void* looper(struct ThreadData* data) {

unsigned int start = data->start;
unsigned int stop = data->stop;
unsigned long N = data->N;
double *R_pairs = data->R_pairs;
double *v = data->v;
unsigned int N_eigs = data->N_eigs;
int scores_as_entries = data->scores_as_entries;
double *pairs_scores = data->pairs_scores;
double *v2 = data->v2;
int (*ALTS)[4][3] = data->ALTS;
int (*signs_confs)[3] = data->signs_confs;

unsigned long N_pairs = N*(N-1)/2;
unsigned int i,j,k,ij,jk,ik;
double s_ij_jk,s_ij_ik,s_ik_jk;
double *Rij, *Rjk, *Rik;
double JRijJ[9], JRjkJ[9], JRikJ[9];
double c[4];
double tmp[9];
int best_i;
double best_val, alt_ij_jk, alt_ik_jk, alt_ij_ik;
unsigned int eig;
double triangle_pairs_score;

for (i=start; i<stop; i++) {
	for (j=i+1; j<(N-1); j++) {
		ij = PAIR_IDX(N,i,j);

		for (k=j+1; k<N; k++) {
			jk = PAIR_IDX(N,j,k);
			ik = PAIR_IDX(N,i,k);
			
			/* compute configurations matches scores */
			Rij = R_pairs + 9*ij;
			Rjk = R_pairs + 9*jk;
			Rik = R_pairs + 9*ik;
			
			JRJ(Rij, JRijJ);
			JRJ(Rjk, JRjkJ);
			JRJ(Rik, JRikJ);
			
			mult_3x3(tmp,Rij,Rjk);
			c[0] = diff_norm_3x3(tmp,Rik);
			
			mult_3x3(tmp,JRijJ,Rjk);
			c[1] = diff_norm_3x3(tmp,Rik);
			
			mult_3x3(tmp,Rij,JRjkJ);
			c[2] = diff_norm_3x3(tmp,Rik);
			
			mult_3x3(tmp,Rij,Rjk);
			c[3] = diff_norm_3x3(tmp,JRikJ);
			
			/* find best match */
			best_i=0; best_val=c[0];
			if (c[1]<best_val) {best_i=1; best_val=c[1];}
			if (c[2]<best_val) {best_i=2; best_val=c[2];}
			if (c[3]<best_val) {best_i=3; best_val=c[3];}
			
			if (scores_as_entries == 1) {
                /* for each triangle side, find the best alternative */
                alt_ij_jk = c[ALTS[0][best_i][0]];
                if (c[ALTS[1][best_i][0]] < alt_ij_jk) {alt_ij_jk = c[ALTS[1][best_i][0]];}
                alt_ik_jk = c[ALTS[0][best_i][1]];
                if (c[ALTS[1][best_i][1]] < alt_ik_jk) {alt_ik_jk = c[ALTS[1][best_i][1]];}
                alt_ij_ik = c[ALTS[0][best_i][2]];
                if (c[ALTS[1][best_i][2]] < alt_ij_ik) {alt_ij_ik = c[ALTS[1][best_i][2]];}

                /* compute triangles entries to be triangles scores */
                s_ij_jk = signs_confs[best_i][0] * (1 - sqrt(best_val/alt_ij_jk));
                s_ik_jk = signs_confs[best_i][1] * (1 - sqrt(best_val/alt_ik_jk));
                s_ij_ik = signs_confs[best_i][2] * (1 - sqrt(best_val/alt_ij_ik));
            }
            else {
                if (scores_as_entries == 2) {
                    /* set triangles entries according to pairs scores */
                    triangle_pairs_score = pairs_scores[ij]*pairs_scores[jk]*pairs_scores[ik];
                    s_ij_jk = signs_confs[best_i][0] * triangle_pairs_score;
                    s_ik_jk = signs_confs[best_i][1] * triangle_pairs_score;
                    s_ij_ik = signs_confs[best_i][2] * triangle_pairs_score;
                }
                else { /* scores_as_entries == 0 */
                    /* set triangles entries to be signs */
                    s_ij_jk = signs_confs[best_i][0];
                    s_ik_jk = signs_confs[best_i][1];
                    s_ij_ik = signs_confs[best_i][2];
                }
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

return NULL;

}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
/* main */

/* output */
#define V_OUT plhs[0]
#define CORES_USED plhs[1]

/* input */
unsigned long N = mxGetScalar(prhs[0]);
double *R_pairs = mxGetPr(mxDuplicateArray(prhs[1]));
double *v = mxGetPr(mxDuplicateArray(prhs[2]));
unsigned int N_eigs = mxGetScalar(prhs[3]);
int scores_as_entries = mxGetScalar(prhs[4]);
double *pairs_scores = mxGetPr(mxDuplicateArray(prhs[5]));
unsigned int n_threads = mxGetScalar(prhs[6]);

/* initialization */
unsigned long N_pairs = N*(N-1)/2;
unsigned long long N_triplets = N*(N-1)*(N-2)/6;
V_OUT = mxCreateDoubleMatrix(N_eigs*N_pairs,1,mxREAL);
double *v2 = mxGetPr(V_OUT);
CORES_USED = mxCreateDoubleMatrix(1,1,mxREAL);
double *cores_used = mxGetPr(CORES_USED);
double **split_v2;
unsigned int i,k,start,stop;
int signs_confs[4][3], ALTS[2][4][3];
unsigned long long tasks_per_thread, n_tasks;
unsigned int n_threads_max;
struct ThreadData *data;
pthread_t *thread;

/* initialize signs configurations and alternatives */
for(i=0; i<4; i++) { for(k=0; k<3; k++) { signs_confs[i][k]=1; } }
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

/* choose number of threads */
n_threads_max = sysconf(_SC_NPROCESSORS_ONLN);
if (n_threads <= 0 || n_threads > n_threads_max) {
	n_threads = n_threads_max;
}
cores_used[0] = (double) n_threads;
tasks_per_thread = (N_triplets+n_threads-1)/n_threads;

/* allocate arrays */
data = (struct ThreadData *) mxMalloc(sizeof(struct ThreadData)*n_threads);
thread = (pthread_t *) mxMalloc(sizeof(pthread_t)*n_threads);
split_v2 = (double **) mxMalloc(sizeof(double *)*n_threads);

/* initialize ThreadData */
stop = 0;
for (k=0; k<n_threads; k++) {
	/* copy basic data */
	data[k].N = N;
	data[k].R_pairs = R_pairs;
	data[k].v = v;
	data[k].N_eigs = N_eigs;
	data[k].scores_as_entries = scores_as_entries;
	data[k].pairs_scores = pairs_scores;
	data[k].ALTS = ALTS;
	data[k].signs_confs = signs_confs;
	
	/* allocate split output arrays */
	split_v2[k] = (double *) mxMalloc(sizeof(double)*N_eigs*N_pairs);
	for (i=0; i<N_eigs*N_pairs; i++) {split_v2[k][i] = 0;}
	data[k].v2 = split_v2[k];
	
	/* compute loop limits */
	start = stop;
	n_tasks = 0;
	for (i=start; i<N-2; i++) {
		n_tasks += (N-i-1)*(N-i-2)/2;
		if (n_tasks>=tasks_per_thread) {break;}
	}
	i++;
	stop = i; if (stop>N-2) {stop=N-2;}
	data[k].start = start;
	data[k].stop = stop;
}

/* launch threads */
for (k=0; k<n_threads; k++) {
	pthread_create(&thread[k], NULL, looper, &data[k]);
}

for (k=0; k<n_threads; k++) {
	pthread_join(thread[k], NULL);
}

/* summarize results */
for (k=0; k<n_threads; k++) {
	for (i=0; i<N_eigs*N_pairs; i++) {
		v2[i] += split_v2[k][i];
	}
}

return;

}
