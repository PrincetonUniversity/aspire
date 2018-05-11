/*

Pairs Probabilities
[ln_f_ind, ln_f_arb, cores_used] = pairs_probabilities_mex...
    (N, R_pairs, P2,A,a,B,b,x0, n_threads, verbose);

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
The variable
	(m_2nd_best - m_best) / m_2nd_best
which represents for every triangle the significance of the best configuration's match,
has an empirically-known approximated distribution, which is divided into two
distributions - one of "good" triangles and one of "bad" triangles.
According to the values of this variable for all the triangles {ijk} of a pair ij,
one can estimate the probability of the pair to have an indicative common line,
through a Bayesian probabilistic model.
This function computes the probability of a pair ij to have the observed values of
these triangles scores - once given it has an indicative common line, and once given
it has an arbitrary common line. The probability of the common line to be indicative
can be derived by Bayes Theorem.

Input:
N - number of projections
R_pairs - relative rotations matrices (3 X 3 X N-choose-2)
P2,A,a,B,b,x0 - scalar parameters of the probabilistic model
n_threads - how many cores to use (if available)
verbose - how detailed to print information - unused by now

Output:
ln_f_ind - ln(the probability of a pair ij to have the observed triangles scores,
given it has an indicative common line) (N-choose-2 X 1)
ln_f_arb - ln(the probability of a pair ij to have the observed triangles scores,
given it has an arbitrary common line) (N-choose-2 X 1)
cores_used - how many cores were used in practice

Compilation (with inline functions):
mex CFLAGS="\$CFLAGS -std=c99" pairs_probabilities_mex.c

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
	double P2,A,a,B,b,x0;
	double *ln_f_ind;
	double *ln_f_arb;
	int (*ALTS)[4][3];
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
double P2 = data->P2;
double A = data->A;
double a = data->a;
double B = data->B;
double b = data->b;
double x0 = data->x0;
double *ln_f_ind = data->ln_f_ind;
double *ln_f_arb = data->ln_f_arb;
int (*ALTS)[4][3] = data->ALTS;

unsigned int i,j,k,ij,jk,ik;
double s_ij_jk,s_ij_ik,s_ik_jk;
double *Rij, *Rjk, *Rik;
double JRijJ[9], JRjkJ[9], JRikJ[9];
double c[4];
double tmp[9];
int best_i;
double best_val, alt_ij_jk, alt_ik_jk, alt_ij_ik;
double f_ij_jk, f_ik_jk, f_ij_ik;

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
			
			/* update probabilities calculations */
            
            /* the probability of a pair ij to have the observed triangles scores,
            given it has an indicative common line */
			f_ij_jk = log( P2*(B*pow(1-s_ij_jk,b)*exp(-b/(1-x0)*(1-s_ij_jk))) + (1-P2)*A*pow((1-s_ij_jk),a) );
			f_ik_jk = log( P2*(B*pow(1-s_ik_jk,b)*exp(-b/(1-x0)*(1-s_ik_jk))) + (1-P2)*A*pow((1-s_ik_jk),a) );
			f_ij_ik = log( P2*(B*pow(1-s_ij_ik,b)*exp(-b/(1-x0)*(1-s_ij_ik))) + (1-P2)*A*pow((1-s_ij_ik),a) );
			ln_f_ind[ij] += f_ij_jk + f_ij_ik;
			ln_f_ind[jk] += f_ij_jk + f_ik_jk;
			ln_f_ind[ik] += f_ik_jk + f_ij_ik;
			
			/* the probability of a pair ij to have the observed triangles scores,
            given it has an arbitrary common line */
			f_ij_jk = log( A*pow((1-s_ij_jk),a) );
			f_ik_jk = log( A*pow((1-s_ik_jk),a) );
			f_ij_ik = log( A*pow((1-s_ij_ik),a) );
			ln_f_arb[ij] += f_ij_jk + f_ij_ik;
			ln_f_arb[jk] += f_ij_jk + f_ik_jk;
			ln_f_arb[ik] += f_ik_jk + f_ij_ik;
			
		}
	}
}

return NULL;

}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
/* main */

/* output */
#define LN_F_IND plhs[0]
#define LN_F_ARB plhs[1]
#define CORES_USED plhs[2]

/* input */
unsigned long N = mxGetScalar(prhs[0]);
double *R_pairs = mxGetPr(mxDuplicateArray(prhs[1]));
double P2 = mxGetScalar(prhs[2]);
double A = mxGetScalar(prhs[3]);
double a = mxGetScalar(prhs[4]);
double B = mxGetScalar(prhs[5]);
double b = mxGetScalar(prhs[6]);
double x0 = mxGetScalar(prhs[7]);
unsigned int n_threads = mxGetScalar(prhs[8]);
int verbose = mxGetScalar(prhs[9]);

/* initialization */
unsigned long N_pairs = N*(N-1)/2;
unsigned long long N_triplets = N*(N-1)*(N-2)/6;
LN_F_IND = mxCreateDoubleMatrix(N_pairs,1,mxREAL);
double *ln_f_ind = mxGetPr(LN_F_IND);
LN_F_ARB = mxCreateDoubleMatrix(N_pairs,1,mxREAL);
double *ln_f_arb = mxGetPr(LN_F_ARB);
CORES_USED = mxCreateDoubleMatrix(1,1,mxREAL);
double *cores_used = mxGetPr(CORES_USED);
double **split_ln_f_ind;
double **split_ln_f_arb;
unsigned int i,k,start,stop;
int ALTS[2][4][3];
unsigned long long tasks_per_thread, n_tasks;
unsigned int n_threads_max;
struct ThreadData *data;
pthread_t *thread;

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
for (i=0; i<N_pairs; i++) {ln_f_ind[i] = 0;}
for (i=0; i<N_pairs; i++) {ln_f_arb[i] = 0;}

/* choose number of threads */
n_threads_max = sysconf(_SC_NPROCESSORS_ONLN);
if (n_threads <= 0 || n_threads > n_threads_max) {
	n_threads = n_threads_max;
}
cores_used[0] = (double) n_threads;
tasks_per_thread = (N_triplets+n_threads-1)/n_threads;
if (verbose >= 1) {
	mexPrintf("Number of available cores: %d\n", n_threads_max);
	mexPrintf("Number of cores used: %d\n", n_threads);
}

/* allocate arrays */
data = (struct ThreadData *) mxMalloc(sizeof(struct ThreadData)*n_threads);
thread = (pthread_t *) mxMalloc(sizeof(pthread_t)*n_threads);
split_ln_f_ind = (double **) mxMalloc(sizeof(double *)*n_threads);
split_ln_f_arb = (double **) mxMalloc(sizeof(double *)*n_threads);

/* initialize ThreadData */
stop = 0;
for (k=0; k<n_threads; k++) {
	/* copy basic data */
	data[k].N = N;
	data[k].R_pairs = R_pairs;
	data[k].P2 = P2;
	data[k].A = A;
	data[k].a = a;
	data[k].B = B;
	data[k].b = b;
	data[k].x0 = x0;
	data[k].ALTS = ALTS;
	
	/* allocate split output arrays */
	split_ln_f_ind[k] = (double *) mxMalloc(sizeof(double)*N_pairs);
	for (i=0; i<N_pairs; i++) {split_ln_f_ind[k][i] = 0;}
	data[k].ln_f_ind = split_ln_f_ind[k];
	split_ln_f_arb[k] = (double *) mxMalloc(sizeof(double)*N_pairs);
	for (i=0; i<N_pairs; i++) {split_ln_f_arb[k][i] = 0;}
	data[k].ln_f_arb = split_ln_f_arb[k];
	
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
	for (i=0; i<N_pairs; i++) {
		ln_f_ind[i] += split_ln_f_ind[k][i];
		ln_f_arb[i] += split_ln_f_arb[k][i];
	}
}

return;

}
