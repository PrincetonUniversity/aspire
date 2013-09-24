// external
extern void dgemm_();
extern void dgemv_();
extern void dgeqp3_();
extern void dsymv_();
extern void dsymm_();
extern void dsyr_();
extern void dsyrk_();
extern void dsyr2k_();
extern void dpotrf_();
extern void dtrsv_();
extern void dtrsm_();
extern void dsaupd_();
extern void dseupd_();
extern void dsyev_();
extern void dtrmm_();
extern size_t  idamax_();
extern size_t gsl_poly_solve_cubic(double, double, double, double*, double*, double*);
extern double gsl_poly_eval(double*, size_t, double);
extern size_t    daxpy_();
extern size_t    dcopy_();
extern double ddot_();
extern double dnrm2_();
extern size_t    dscal_();

/* extern 
        dgemm_(&transr, &transv, &rk, &(A->lr->ncol), &(A->lr->nrow), &one,
               U + base + 1, &(A->lr->nrow), A->lr->ent + 1, &(A->lr->nrow), &zero,
               global_UtB + 1, &rk); */
/*extern void dgemm_(char *, char *, size_t *, size_t *, size_t *, double *, double *, size_t *, double *b, size_t *, double *, double *, size_t *ldc);*/
