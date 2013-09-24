/* copystructures.c */
size_t copystructures(problemdata *data, size_t m, size_t numblk, size_t *blksz, char *blktype, double *b, double *CAent, size_t *CArow, size_t *CAcol, size_t *CAinfo_entptr, char *CAinfo_type);
/* dataoper.c */
double function(problemdata *data, double *R);
size_t Aoper(problemdata *d, double *U, double *V, double *UVt, size_t same, size_t obj, double *results);
size_t Aoper_formUVt(problemdata *data, double *passedUVt, double *U, double *V, size_t same);
size_t gradient(problemdata *data, double *R);
size_t StimesR(problemdata *data, double *S, double *y, double *R, double *result);
size_t Stimesmat(problemdata *data, double *S, double *y, double *vec, double *result, size_t n, size_t m, size_t k);
size_t AToper(problemdata *data, double *y, double *S, size_t obj);
double C_normdatamat(problemdata *data);
double normdatamat(problemdata *data, size_t matnum);
size_t essential_calcs(problemdata *data, double *R, double normC, double normb, double *val, double *rho_c_val, double *rho_f_val);
/* eigen.c */
size_t my_arpack(size_t (*matvec)(double *, double *), size_t n, size_t nev, size_t ncv, char *which, size_t maxitr, size_t printlevel, double *evals, double *evecs, size_t *nconv, size_t *nummatvec);
size_t simple_Stimesvec_block(double *out, double *in);
int Smineval(problemdata *d, double *eval);
int Sblockmineval(problemdata *d, double *blkevals);
/* initialize.c */
size_t initialize(problemdata *data, size_t *maxranks);
size_t deinitialize(problemdata *data);
size_t locatetype(problemdata *data, size_t blk, char type, size_t **passed_ind, size_t *passed_nind);
/* lbfgs.c */
size_t dirlbfgs(problemdata *data, lbfgsvec *vecs, double *dir, double *grad, size_t oldest, size_t numlbfgsvecs, size_t scale);
size_t updatelbfgs1(problemdata *data, lbfgsvec *vecs, double *grad, size_t oldest);
size_t updatelbfgs2(problemdata *data, lbfgsvec *vecs, double *dir, double *grad, double stepsize, size_t *oldest, size_t numlbfgsvecs);
/* linesearch.c */
double linesearch(problemdata *data, double *R, double *D, double max, double *funcval, size_t update);
/* main.c */
int main(size_t argc, char *argv[]);
size_t getstorage(size_t m, size_t numblk, size_t *blksz, char *blktype, size_t *CAinfo_entptr, size_t *passedn, size_t *passednr, size_t *maxranks);
size_t readin(size_t m, size_t numblk, size_t *blksz, char *blktype, double *R, double *lambda, size_t *maxranks, size_t *ranks, double *pieces, FILE *fid);
size_t writeout(size_t m, size_t numblk, size_t *blksz, char *blktype, double *R, double *lambda, size_t *maxranks, size_t *ranks, double *pieces, FILE *fid);
/* mexsdplr.c */
size_t classifyApq(size_t *blksz, char *blktype, size_t *inblk, size_t *inblkptr, size_t p, size_t q, size_t A_tr, size_t *cons, size_t *blk, size_t *i, size_t *j);
/* misc.c */
size_t print_notes(void);
size_t printparams(problemdata *data);
size_t print_dimacs_errors(problemdata *data, double *R);
size_t printheading(size_t start);
/* params.c */
size_t getparams(char *paramfile, size_t *inputtype, double *rho_f, double *rho_c, double *sigmafac, size_t *rankreduce, size_t *timelim, size_t *printlevel, size_t *dthresh_dim, double *dthresh_dens, size_t *numbfgsvecs, double *rankredtol, double *gaptol, ptrdiff_t *checkbd, size_t *typebd);
size_t getparams_maxlinelength(FILE *datafile);
size_t getparams_getline(FILE *datafile, char *buffer, size_t bufsiz);
size_t getparams_tolower(char *buff, size_t buffsz);
size_t getparams_isnumber(char *str);
size_t generate_params(void);
size_t generate_params_info(size_t num);
/* rankreduce.c */
size_t dorankreduce(problemdata *d, double *R);
/* readdata.c */
size_t readdata_sdpa(char *datafilename, size_t *passed_m, size_t *passed_numblk, ptrdiff_t **passed_blksz, char **passed_blktype, double **passed_b, double **passed_CAent, size_t **passed_CArow, size_t **passed_CAcol, size_t **passed_CAinfo_entptr, size_t **passed_CAinfo_rowcolptr, char **passed_CAinfo_type, char **passed_CAinfo_storage);
size_t quicksort5(size_t *A1, size_t *A2, size_t *A3, size_t *A4, double *A5, size_t p, size_t r);
size_t partition5(size_t *A1, size_t *A2, size_t *A3, size_t *A4, double *A5, size_t p, size_t r);
void skip_to_end_of_line(FILE *datafile);
size_t get_line(FILE *datafile, char *buffer, size_t bufsiz);
size_t max_line_length(FILE *datafile);
size_t readdata_sdplr(char *datafilename, size_t *passed_m, size_t *passed_numblk, ptrdiff_t **passed_blksz, char **passed_blktype, double **passed_b, double **passed_CAent, size_t **passed_CArow, size_t **passed_CAcol, size_t **passed_CAinfo_entptr, size_t **passed_CAinfo_rowcolptr, char **passed_CAinfo_type, char **passed_CAinfo_storage);
size_t readdata_raw(char *datafilename, size_t *passed_m, size_t *passed_numblk, ptrdiff_t **passed_blksz, char **passed_blktype, double **passed_b, double **passed_CAent, size_t **passed_CArow, size_t **passed_CAcol, size_t **passed_CAinfo_entptr, size_t **passed_CAinfo_rowcolptr, char **passed_CAinfo_type, char **passed_CAinfo_storage);
size_t writedata_raw(char *datafilename, size_t m, size_t numblk, ptrdiff_t *blksz, char *blktype, double *b, double *CAent, size_t *CArow, size_t *CAcol, size_t *CAinfo_entptr, size_t *CAinfo_rowcolptr, char *CAinfo_type, char *CAinfo_storage);
size_t writedata_sdplr(char *datafilename, size_t m, size_t numblk, ptrdiff_t *blksz, char *blktype, double *b, double *CAent, size_t *CArow, size_t *CAcol, size_t *CAinfo_entptr, char *CAinfo_type);
size_t writedata_sdpa(char *datafilename, size_t m, size_t numblk, ptrdiff_t *blksz, char *blktype, double *b, double *CAent, size_t *CArow, size_t *CAcol, size_t *CAinfo_entptr, char *CAinfo_type);
/* sdplrlib.c */
size_t sdplrlib(size_t m, size_t numblk, ptrdiff_t *blksz, char *blktype, double *b, double *CAent, size_t *CArow, size_t *CAcol, size_t *CAinfo_entptr, char *CAinfo_type, size_t numbfgsvecs, double rho_f, double rho_c, double sigmafac, size_t rankreduce, double gaptol, ptrdiff_t checkbd, size_t typebd, size_t dthresh_dim, double dthresh_dens, size_t timelim, double rankredtol, size_t printlevel, double *R, double *lambda, size_t *maxranks, size_t *ranks, double *pieces);
void myprint(size_t majiter, size_t iter, double val, double rho_f_val, size_t CG, double totaltime);
size_t do_scaling(problemdata *data, double value, double *norm);
/* timefuncs.c */
#ifndef __WIN32
double current_time(double timeorig);
#else
double current_time(clock_t timeorig);
#endif
/* util.c */
size_t copyscaledvectovec(double *dy, double da, double *dx, size_t n);
size_t move_in_dir(double *dy, double *dx, double da, double *dir, size_t n);
size_t mydaxpy(size_t n, double da, double *dx, size_t incx, double *dy, size_t incy);
size_t mydcopy(size_t n, double *dx, size_t incx, double *dy, size_t incy);
double myddot(size_t n, double *dx, size_t incx, double *dy, size_t incy);
double mydnrm2(size_t n, double *dx, size_t incx);
size_t mydscal(size_t n, double da, double *dx, size_t incx);
size_t createlowrankmat(lowrankmat **passedR, size_t ncol, size_t nrow);
size_t destroylowrankmat(lowrankmat *R);
size_t createsparsesymmmat(sparsesymmmat **passedS, size_t nnz);
size_t destroysparsesymmmat(sparsesymmmat *S);
size_t creatediagmat(diagmat **passedD, size_t nnz);
size_t destroydiagmat(diagmat *D);
size_t createdatamat(datamat **passedA, char type, size_t ncol_or_nnz, size_t dim, char *label);
size_t destroydatamat(datamat *A);
