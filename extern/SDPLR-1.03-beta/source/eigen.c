#include "myinclude.h"
#ifdef __MEX
#include "mex.h"
#endif

/* #ifdef __MEX */
/* #undef exit */
/* #define exit(er) mexWarnMsgTxt("SDPLR Error\n"); return mxGetNaN(); */
/* #endif */

#ifdef __ARPACK

size_t my_arpack(size_t (*matvec)(double*, double*), size_t n, size_t nev, size_t ncv, char* which, size_t maxitr, size_t printlevel, double* evals, double* evecs, size_t* nconv, size_t* nummatvec)
{
  // Variables taken from ARPACK examples
  double *v, *workl, *workd, *d, *resid;
  size_t *select, iparam[11], ipntr[11];
  char bmat[2], all[4];
  size_t ldv, ido, lworkl, info, ierr;
  size_t mode, ishfts, rvec;
  double tol, sigma;
  double zero;

  size_t i;
  double minevalest;

  // Extra variables
  size_t tempsize_t;

  // Initialize counter for number of matvec operations
  *nummatvec = 0;

  // Setup different parameters
  ldv = n;
  zero = 0.0;
  strcpy(bmat, "I\0");
  lworkl = ncv*(ncv+8);
  //tol = zero;
  tol = (double)1.0/n;
  info = 0;
  ido = 0;
  ishfts = 1;
  mode = 1;
  iparam[1-1] = ishfts;
  iparam[3-1] = maxitr;
  iparam[7-1] = mode;
  rvec = 1;
  strcpy(all, "All\0");

  // Allocate memory
  v = (double*)calloc(ldv*ncv, sizeof(double));
  workl = (double*)calloc(ncv*(ncv+8), sizeof(double));
  workd = (double*)calloc(3*n, sizeof(double));
  d = (double*)calloc(ncv*2, sizeof(double));
  resid = (double*)calloc(n, sizeof(double));
  select = (size_t*)calloc(ncv, sizeof(size_t));

  // Do ARPACK loop

  do {

    dsaupd_(&ido, bmat, &n, which, &nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, &info);

    if(ido == -1 || ido == 1) {

      // This line assumes matvec numbers from 0 (not 1)
      matvec(workd + ipntr[2-1] - 1, workd + ipntr[1-1] - 1);
      *nummatvec = *nummatvec + 1;

    }
    else if(info < 0) {

      if(printlevel > 0) printf("ARPACK error (info = %d).\n", info);
      goto END_ROUTINE;

    }

  } while(ido != 99);

  // Post process

  dseupd_(&rvec, all, select, d, v, &ldv, &sigma, bmat, &n, which, &nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, &ierr);

  if(ierr != 0) {
    if(printlevel > 0) printf("ARPACK error (ierr = %d).\n", ierr);
    info = ierr;
    goto END_ROUTINE;
  }

  *nconv = iparam[5-1];

  if(printlevel > 0 && *nconv < nev) {
    printf("Warning: %d out of %d evals converged.", *nconv, nev);
    if(nev == 1) printf(" Best estimate returned.\n");
    else printf("\n");
  }

  // Save eigen pairs in user's memory (even unreliable ones)
  mydcopy(nev, d, 1, evals, 1);
  tempsize_t = nev*n;
  mydcopy(tempsize_t, v, 1, evecs, 1);

  if(nev == 1 && *nconv == 0) {
    minevalest = 1.0e10;
    for(i = 1; i <= ncv; i++)
      minevalest = mymin(minevalest, workl[ipntr[6-1]-2+i]);
    evals[0] = minevalest;
  }


END_ROUTINE:

  free(select);
  free(resid);
  free(d);
  free(workd);
  free(workl);
  free(v);


  if(*nconv < nev) return -1;
  return info;
}

#endif

size_t simple_Stimesvec_block(double* out, double* in)
{
  // This routine is called by ARPACK, which assumes
  // out and in are numbered from 0

  Stimesmat(global_data, global_data->S + global_data->XS_blkptr[global_blk] - 1, global_data->y, in - 1, out - 1, global_data->blksz[global_blk], 1, global_blk);

  return 0;
}

int Smineval(problemdata* d, double* eval)
{
  int problem, ct, i, k;
  double *blkevals;

  problem = 0;

  *eval = 1.0e10;

  ct = 0;
  for(k = 1; k <= d->numblk; k++)
    if(d->blktype[k] == SDPBLK) ct++;
    else if(d->blktype[k] == DIAGBLK) ct += d->blksz[k];

  MYCALLOC(blkevals, double, ct+1);

  problem = Sblockmineval(d, blkevals);

  for(i = 1; i <= ct; i++)
    *eval = mymin(*eval, blkevals[i]);

  MYFREE(blkevals);

  return problem;
}

int Sblockmineval(problemdata* d, double* blkevals)
{
  int problem=0;
  size_t ct, i, k, nev, maxitr, printlevel;
  char which[3];
  double eval;
  size_t maxblksz, info, dimwork;
  char uplo, jobz;
  double *work, *Scopy, *evals;

  // Note this function overwrites S!

  
  global_data = d;
  nev = 1;
  strcpy(which, "SA\0");
  maxitr = 500;
  printlevel = 0;


  // If necessary, allocate workspace for dsyev_
  maxblksz = 0;
  for(k = 1; k <= d->numblk; k++)
    if(d->blktype[k] == SDPBLK && d->XS_blksto[k] == DENSE)
      maxblksz = mymaxint(maxblksz, d->blksz[k]);
  if(maxblksz > 0) {
    dimwork = mymaxint(1, 3*maxblksz - 1);
    MYCALLOC(Scopy, double, maxblksz*maxblksz + 1);
    MYCALLOC(work, double, dimwork + 1);
    MYCALLOC(evals, double, maxblksz + 1);
  }
  else { work = Scopy = evals = NULL; }

  ct = 0;

  for(k = 1; k <= d->numblk; k++) {

    if(d->blktype[k] == SDPBLK && d->XS_blksto[k] == SPARSE) {
      
#ifdef __ARPACK

      global_blk = k;
      n = d->blksz[k];
      ncv = mymin(2*mymax(2*nev,20),n); // This is Matlab's rule (times two).
      MYCALLOC(evec, double, n+1);

      // When finished, evals[1] gives first eigenvalue
      // and evecs[1] is start of eigenvectors
      ret = my_arpack(simple_Stimesvec_block, n, nev, ncv, which, maxitr, printlevel, &eval, evec+1, &nconv, &nummatvec);
      MYFREE(evec);
      if(ret == -1) problem = -1;

#else

      eval = -1.0e10;

#endif

      blkevals[++ct] = eval;

    }
    else if(d->blktype[k] == SDPBLK && d->XS_blksto[k] == DENSE) {

        jobz = 'n'; uplo = 'l';
        EASYDCOPY(d->blksz[k]*d->blksz[k], d->S + d->XS_blkptr[k] - 1, Scopy);
        dsyev_(&jobz, &uplo, &(d->blksz[k]), Scopy + 1, &(d->blksz[k]), evals + 1, work + 1, &dimwork, &info);
        if((int)info != 0) { printf("Eigenvalue computation failed.\n"); exit(0); }
        blkevals[++ct] = evals[1]; 

    }
    else if(d->blktype[k] == DIAGBLK) {

      for(i = 1; i <= d->blksz[k]; i++)
        blkevals[++ct] = d->S[d->XS_blkptr[k] - 1 + i];

    }

  }

  // If necessary, deallocate workspace used for dsyev_
  if(maxblksz > 0) {
    MYFREE(Scopy);
    MYFREE(work);
    MYFREE(evals);
  }

  return problem;
}

// double dualbound(problemdata* data)
// {
//   size_t problem;
//   double eval;

//   EASYDCOPY(data->m, data->lambda, data->S[1]->y);
//   problem = Smineval(data, &eval);
//  
//   return EASYDDOT(data->m, data->b, data->S[1]->y) - data->evaladj*mymin(eval, 0.0);

// }


