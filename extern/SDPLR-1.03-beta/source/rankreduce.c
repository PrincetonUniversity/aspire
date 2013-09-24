#include "myinclude.h"
#ifdef __MEX
#include "mex.h"
#endif

size_t dorankreduce(problemdata* d, double* R)
{
  size_t h, i, j, k, p, q, base, base1, base2;
  size_t maxrank, maxblksz, newrank=0, **save;
  double tol, *U, *tempU, tv;
  size_t info, lwork, *jpvt;
  double *work, *tau;
  double currnorm, totalnorm;

  tol = d->rankredtol; // machine precision (double) for P4 as given by Matlab

  // Calculate maxrank and maxblksz (for storage purposes)
  maxrank = -1; maxblksz = -1;
  for(k = 1; k <= d->numblk; k++)
    if(d->blktype[k] == 's') {
      maxrank = mymax(maxrank, d->rank[k]);
      maxblksz = mymax(maxblksz, d->blksz[k]);
    }

  // Allocate space that will be need in this subroutine
  MYCALLOC(U,     double, d->nr     + 1);
  MYCALLOC(tempU, double, d->nr     + 1);
  MYCALLOC(save,  size_t*,   d->numblk + 1);
  for(k = 1; k <= d->numblk; k++)
    MYCALLOC(save[k], size_t, d->rank[k] + 1);

  lwork = mymax(1, 3*maxblksz + 1);
  MYCALLOC(work, double, lwork    + 1);
  MYCALLOC(jpvt, size_t,    maxblksz + 1);
  MYCALLOC(tau,  double, maxrank  + 1);


  // Go through SDP blocks to get numerical ranks,
  // as well as basis for range space of R. Determine
  // which "columns" of R will be saved.

  base = 0;

  for(k = 1; k <= d->numblk; k++) {

    if(d->blktype[k] == 's') {

      p = d->rank[k];
      q = d->blksz[k];

      for(i = 1; i <= q; i++) jpvt[i] = 0;

      // copy R^T size_to U
      for(i = 1; i <= q; i++) for(j = 1; j <= p; j++)
          U[base + (i-1)*p + j] = R[base + (j-1)*q + i];

      // Do QR factorization of U
      dgeqp3_(&p, &q, U + base + 1, &p, jpvt + 1, tau + 1, work + 1, &lwork, &info);
      if(info != 0) {
        printf("Error (dorankreduce.c): Problem with QR factorization! (info = %d)\n", info);
        exit(0);
      }

      // Overwrite strictly lower part of U with 0 (this is necessary later)
      for(i = 1; i <= p; i++) for(j = 1; j <= i-1; j++)
          U[base + (j-1)*p + i] = 0.0;

      // Go through diagonal of U, and determine which columns fit the tolerance (save[k][i] = 1)
      totalnorm = EASYDNRM2(p*q, U + base);
      currnorm = 0.0;
// OLD
//       for(i = p; i >= 1; i--) {
//         currnorm += sqrt(pow(currnorm, 2.0) + pow(mydnrm2(q, U + base + i, p), 2.0));
//         if(currnorm <= tol*totalnorm) save[k][i] = 0;
//         else save[k][i] = 1;
//       }
      for(i = p; i >= 1; i--) {
        currnorm += pow(mydnrm2(q, U + base + i, p), 2.0);
        if(currnorm <= (tol*tol)*(totalnorm*totalnorm)) save[k][i] = 0;
        else save[k][i] = 1;
      }
      currnorm = sqrt(currnorm);

      // Double check that cols which fit tolerance are at beginning of list;
      // this is a feature of dgeqp3_().
      // Apparent problems b/c of roundoff. For safety, will just force cols
      // that we are unsure of to be saved.
      info = 0;
      for(i = 2; i <= p; i++) if(save[k][i] != save[k][i-1]) info++;
      if(info >= 2) {
        //info = 0;
        //for(i = p; i >= 1; i--) {
        //  if(save[k][i] == 1) info = 1;
        //  save[k][i] = info;
        //}
        printf("Error (rankreduce.c): Unexpected results from QR factorization.\n");
        //for(i = 1; i <= p; i++)
        //  printf("%e\n", U[base + (i-1)*p + i]);
        exit(0);
      }

      // Copy permuted U size_to tempU and do transpose at the same time
      for(j = 1; j <= q; j++)
        for(i = 1; i <= p; i++)
          tempU[base + (i-1)*q + jpvt[j]] = U[base + (j-1)*p + i];

      // Copy back (then can forget about tempU, and also U is the same orientation as R)
      for(j = 1; j <= q; j++)
        for(i = 1; i <= p; i++)
          U[base + (j-1)*p + i] = tempU[base + (j-1)*p + i];

      // Determine what actually to do (reduce rank or increase rank or nothing);
      // saved in save[k][0] as -1 (reduce) or +1 (increase) or 0 (nothing)
      info = 1;
      for(i = 1; i <= p; i++) info *= save[k][i];
      if(info == 1 && p == d->maxrank[k]) save[k][0] = 0;
      if(info == 1 && p < d->maxrank[k]) save[k][0] = 1;
      if(info == 0 && p == 1) save[k][0] = 0;
      if(info == 0 && p >= 2 && save[k][p-1] == 0) save[k][0] = -1;
      if(info == 0 && p >= 2 && save[k][p-1] == 1) save[k][0] = 0;

      // Adjust save[k] to indicate what will actually be done in
      // the case of reducing the rank
      if(save[k][0] == -1)
        for(i = 1; i <= p; i++)
          if(save[k][i] == 0) {
            save[k][i] = 1;
            break;
          }

    }
    else if(d->blktype[k] == 'd') {

      EASYDCOPY(d->blksz[k]*d->rank[k], R + base, U + base);

    }

    base += d->blksz[k]*d->rank[k];

  }

  // Now change structure of R

  base1 = 0;
  base2 = 0;

  for(k = 1; k <= d->numblk; k++) {

    if(d->blktype[k] == 's') {

      if(save[k][0] != 0) {

        if(save[k][0] == -1) {

          // Copy U size_to R
          h = 0;
          for(i = 1; i <= d->rank[k]; i++)
            if(save[k][i] == 1) {
              h++;
              EASYDCOPY(d->blksz[k], U + base2 + (i-1)*d->blksz[k], R + base1 + (h-1)*d->blksz[k]);
            }
          newrank = h;

        }
        else if(save[k][0] == 1) {

          EASYDCOPY(d->blksz[k]*d->rank[k], U + base2, R + base1);

          if(RANDOM) srand( (unsigned)time( NULL ) );
          else       srand(925);

          tv = 0.0;
          for(h = 1; h <= d->blksz[k]*(d->maxrank[k] - d->rank[k]); h++) {
            R[base1 + d->blksz[k]*d->rank[k] + h]  = (double)rand()/RAND_MAX;
            R[base1 + d->blksz[k]*d->rank[k] + h] -= (double)rand()/RAND_MAX;
            tv += pow(R[base1 + d->blksz[k]*d->rank[k] + h], 2.0);
          }
          tv = sqrt(tv);
          for(h = 1; h <= d->blksz[k]*(d->maxrank[k] - d->rank[k]); h++)
            R[base1 + d->blksz[k]*d->rank[k] + h] /= (double)d->blksz[k]*tv;

          newrank = d->maxrank[k];

        }

      }
      else if(save[k][0] == 0) {

        EASYDCOPY(d->blksz[k]*d->rank[k], U + base2, R + base1);
        newrank = d->rank[k];

      }
        
      //printf("  <rr> blk %d  ( %d / %d / %d )\n", k, newrank, d->rank[k], d->maxrank[k]);

    }
    else if(d->blktype[k] == 'd') {

      EASYDCOPY(d->blksz[k]*d->rank[k], U + base2, R + base1);
      newrank = d->rank[k];

    }

    base2 += d->blksz[k]*d->rank[k];
    d->rank[k] = newrank; // Update rank
    base1 += d->blksz[k]*d->rank[k];

  }

  // Update total dimension nr
  d->nr = base1;


  // Clean up
  MYFREE(U);
  MYFREE(tempU);
  MYFREE(work);
  MYFREE(jpvt);
  MYFREE(tau);
  for(k = 1; k <= d->numblk; k++) MYFREE(save[k]);
  MYFREE(save);

  return 0;

}

