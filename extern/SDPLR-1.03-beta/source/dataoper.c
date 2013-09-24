#include "myinclude.h"
#ifdef __MEX
#include "mex.h"
#endif

/* #ifdef __MEX */
/* #undef exit */
/* #define exit(er) mexWarnMsgTxt("SDPLR Error\n"); return mxGetNaN(); */
/* #endif */

double function(problemdata* data, double* R)
{

  Aoper(data,R,R,data->X,1,1,data->vio);
  EASYDAXPY(data->m, -1.0, data->b, data->vio);

  return data->vio[0] - EASYDDOT(data->m, data->vio, data->lambda)
    + 0.5*data->sigma*pow(EASYDNRM2(data->m, data->vio), 2.0);

}

size_t Aoper(problemdata* d, double* U, double* V, double* UVt, size_t same, size_t obj, double* results)
{
  // This function computes 0.5*A(UV^T + VU^T).
  // UVt will store 0.5*(UV^T + VU^T)
  // same = 1 if U = V
  // obj = 1 if 0.5 * C * (UV^T + VU^T) is to be stored in results[0]

  size_t h, i, j, k;
  size_t base, rk;
  double sum;
  double one=1.0, zero=0.0;
  char transr, transv;
  datamat *A;

  size_t sami, samj, samk;
  FILE *fid;

  Aoper_formUVt(d,UVt,U,V,same);

  // Do those A_i which are diagonal or sparse
  
  for(i = 1-obj; i <= d->m; i++) { // Note 'i = 1-obj'
    results[i] = 0.0;
    for(j = d->AA_rowptr[i]; j <= d->AA_rowptr[i+1]-1; j++)
      results[i] += d->AA_colval_two[j]*UVt[d->AA_colind[j]];
  }

  // Now handle those A_i which are lowrank

  if(same) {

    for(h = 1; h <= d->lr_num; h++) {

      k  = d->lr_blk[h];  
      i  = d->lr_mat[h];  
      rk = d->rank[k];

      base = 0; // Can I fix this?
      for(j = 1; j <= k-1; j++)
        base += d->blksz[j]*d->rank[j];

      if(i > 0 || obj) {

        if(i == 0) A = d->C   [k];
        else       A = d->A[i][k];

        // Assume A_i = B D B^T

        // Form U^T B and V^T B 
        transr = 't'; transv = 'n';
        dgemm_(&transr, &transv, &rk, &(A->lr->ncol), &(A->lr->nrow), &one,
               U + base + 1, &(A->lr->nrow), A->lr->ent + 1, &(A->lr->nrow), &zero,
               global_UtB + 1, &rk);

        // Get value and store it
        sum = 0.0;
        for(j = 1; j <= A->lr->ncol; j++)
          sum += A->lr->d[j]*EASYDDOT(rk, global_UtB + (j-1)*rk, global_UtB + (j-1)*rk);
        results[i] += sum;

      }
    }


  }
  else if(!same) {

    for(h = 1; h <= d->lr_num; h++) {

      k  = d->lr_blk[h];  
      i  = d->lr_mat[h];  
      rk = d->rank[k];

      base = 0; // Can I fix this?
      for(j = 1; j <= k-1; j++)
        base += d->blksz[j]*d->rank[j];

/*      printf("k = %d\n", k);
      printf("i = %d\n", i);
      printf("rk = %d\n", rk);
      printf("base = %d\n", base); */


      if(i > 0 || obj) {

        if(i == 0) A = d->C   [k];
        else       A = d->A[i][k];

/*        fid = fopen("trythis.m","w");
        for(sami = 0; sami < A->lr->nrow; sami++) {
          for(samj = 0; samj < rk; samj++) {
            fprintf(fid, "U(%d,%d) = %.15e;\n", sami+1, samj+1, U[base + 1 + samj*A->lr->nrow + sami]);
          }
        }
        for(sami = 0; sami < A->lr->nrow; sami++) {
          for(samj = 0; samj < rk; samj++) {
            fprintf(fid, "V(%d,%d) = %.15e;\n", sami+1, samj+1, V[base + 1 + samj*A->lr->nrow + sami]);
          }
        }
        printf("norm U = %f\n", EASYDNRM2(A->lr->nrow*rk, U+1));
        printf("norm V = %f\n", EASYDNRM2(A->lr->nrow*rk, V+1)); */


        transr = 't'; transv = 'n';
        dgemm_(&transr, &transv, &rk, &(A->lr->ncol), &(A->lr->nrow), &one,
               U + base + 1, &(A->lr->nrow), A->lr->ent + 1, &(A->lr->nrow), &zero,
               global_UtB + 1, &rk);
        dgemm_(&transr, &transv, &rk, &(A->lr->ncol), &(A->lr->nrow), &one,
               V + base + 1, &(A->lr->nrow), A->lr->ent + 1, &(A->lr->nrow), &zero,
               global_VtB + 1, &rk);

        sum = 0.0;
        for(j = 1; j <= A->lr->ncol; j++)
          sum += A->lr->d[j]*EASYDDOT(rk, global_UtB + (j-1)*rk, global_VtB + (j-1)*rk);
        results[i] += sum;
/*        printf("sum = %f\n", sum);

        exit(0); */

      }
    }

  }

  return 0;

}

//! \page page_Aoper_formUVt Evaluating the Constrasize_ts (Aoper_formUVt in dataoper.c)
//! 
//! This page describes one of the main computational routines in SDPLR,
//! namely to apply the constrasize_t operator \f$ {\cal A} \f$ to matrices
//! of the form \f$ RR^T \f$. In fact, matrices of the more general form
//! \f$ UV^T \f$ are allowed because these are found during the exact
//! \ref linesearch procedure.

size_t Aoper_formUVt(problemdata* data, double* passedUVt, double* U, double* V, size_t same)
{
  size_t i, k, base;
  char uplo='l', trans='n';
  double half=0.5, one=1.0, zero=0.0;

  /* size_t ctr, blksz, rank; */
  size_t blksz, rank;
  /* double tempval; */

  size_t j;
  size_t *colptr, *rowind;
  double *UVt;
  /* double *ptr1, *ptr2, *ptr3, *ptr4; */

  if(same) {

    base = 0;

    for(k = 1; k <= data->numblk; k++) {

      UVt   = passedUVt + data->XS_blkptr[k]-1;
      blksz = data->blksz[k];
      rank  = data->rank[k];

      if(data->blktype[k] == SDPBLK) {

        if(data->XS_blksto[k] == SPARSE) {

          colptr = data->XS_colptr[k];
          rowind = data->XS_rowind[k];

          for(j = 1; j <= blksz; j++)
            for(i = colptr[j]; i <= colptr[j+1]-1; i++)
              UVt[i] = myddot(rank, U + base + rowind[i], blksz, U + base + j, blksz);

        }
        else if(data->XS_blksto[k] == DENSE) {

          dsyrk_(&uplo, &trans, &blksz, &rank, &one, U + base + 1, &blksz, &zero, UVt + 1, &blksz);

        }

      }
      else if(data->blktype[k] == DIAGBLK) {

        for(i = 1; i <= blksz; i++)
          UVt[i] = pow(U[base + i], 2.0);

      }
      else { printf("Aoper_formUVt: Unrecognized blktype.\n"); exit(0); }

      base += blksz*rank;

    }

  }
  else if(!same) {

    base = 0;

    for(k = 1; k <= data->numblk; k++) {

      UVt   = passedUVt + data->XS_blkptr[k]-1;
      blksz = data->blksz[k];
      rank  = data->rank[k];

      if(data->blktype[k] == SDPBLK) {

        if(data->XS_blksto[k] == SPARSE) {

          colptr = data->XS_colptr[k];
          rowind = data->XS_rowind[k];

          for(j = 1; j <= blksz; j++)
            for(i = colptr[j]; i <= colptr[j+1]-1; i++) {
              UVt[i]  = myddot(rank, U + base + rowind[i], blksz, V + base + j, blksz);
              UVt[i] += myddot(rank, V + base + rowind[i], blksz, U + base + j, blksz);
              UVt[i] *= 0.5;
            }

        }
        else if(data->XS_blksto[k] == DENSE) {

          dsyr2k_(&uplo, &trans, &blksz, &rank, &half, U + base + 1, &blksz, V + base + 1, &blksz, &zero, UVt + 1, &blksz);

        }

      }
      else if(data->blktype[k] == DIAGBLK) {

        for(i = 1; i <= blksz; i++)
          UVt[i] = U[base + i]*V[base + i];

      }
      else { printf("Aoper_formUVt: Unrecognized blktype.\n"); exit(0); }

      base += blksz*rank;

    }

  }

  return 1;

}

size_t gradient(problemdata* data, double* R)
{
  size_t i;
  double *G;

  G = data->G;

  data->y[0] = 1.0;
  for(i = 1; i <= data->m; i++)
    data->y[i] = -(data->lambda[i] - data->sigma*data->vio[i]);

  AToper(data, data->y, data->S, 1);

  StimesR(data, data->S, data->y, R, G);

  EASYDSCAL(data->nr, 2.0, G);

  return 1;
}

size_t StimesR(problemdata *data, double *S, double *y, double *R, double *result)
{
  size_t k, base=0;

  for(k = 1; k <= data->numblk; k++) {
    Stimesmat(data, S + data->XS_blkptr[k] - 1, y, R + base, result + base, data->blksz[k], data->rank[k], k);
    base += data->blksz[k]*data->rank[k];
  }

  return 1;
}

size_t Stimesmat(problemdata *data, double *S, double *y, double* vec, double* result, size_t n, size_t m, size_t k)
{
  size_t h, i, j, r;
  size_t *colptr, *rowind;
  datamat *A;
  char side='l', uplo='l', transyes='t', transno='n';
  double one=1.0, zero=0.0, *U;

  // k is the block number

  // both vec and result are n x m

  if(data->blktype[k] == SDPBLK) {

    ZEROVEC(result, n*m);

    // The next part (multiplication by low-rank data matrices) is
    // only done if S is stored as sparse. Otherwise, low-rank data
    // matrices are already incorporated size_to S explicitly.

    if(data->XS_blksto[k] == SPARSE) {

      for(h = 1; h <= data->nlrind[k]; h++) {

        i = data->lrind[k][h];

        if(i == 0) A = data->C   [k];
        else       A = data->A[i][k];

        // Allocate temp space
        MYCALLOC(U, double, A->lr->ncol*m + 1); // Gotta get rid of this

        // Form V^T M, store in U
        dgemm_(&transyes, &transno, &(A->lr->ncol), &m, &(A->lr->nrow), &one,
               A->lr->ent + 1, &(A->lr->nrow), vec + 1, &(A->lr->nrow), &zero,
               U + 1, &(A->lr->ncol));

        // Form DU, store in U
        for(j = 1; j <= A->lr->ncol; j++)
          mydscal(m, A->lr->d[j], U + j, A->lr->ncol);

        // Form VU, put size_to result scaled by ???
        dgemm_(&transno, &transno, &n, &m, &(A->lr->ncol), y+i,
               A->lr->ent + 1, &n, U + 1, &(A->lr->ncol), &one,
               result + 1, &n);

        // Free temp space
        MYFREE(U); // Gotta get rid of this

      }

      // Now do sparse part

      colptr = data->XS_colptr[k];
      rowind = data->XS_rowind[k];

      for(j = 1; j <= data->blksz[k]; j++)
        for(i = colptr[j]; i <= colptr[j+1]-1; i++) {
          r = rowind[i];
                     mydaxpy(m, S[i], vec + r, n, result + j, n);
          if(r != j) mydaxpy(m, S[i], vec + j, n, result + r, n);
        }

    }
    else if(data->XS_blksto[k] == DENSE)
      dsymm_(&side, &uplo, &n, &m, &one, S + 1, &n, vec + 1, &n, &one, result + 1, &n);

  }
  else if(data->blktype[k] == DIAGBLK) {

    for(j = 1; j <= n; j++)
      result[j] = S[j]*vec[j];

  }

  return 1;
  
}

size_t AToper(problemdata* data, double* y, double* S, size_t obj)
{
  size_t      h, i, j, k;
  size_t      inc1 = 1;
  char     uplo = 'l';
  double   origy0=0.0, tv;
  datamat *A;

  // If obj==1, make sure y[0] = 1.0
  
  if(obj) { origy0 = y[0]; y[0] = 1.0; }

  // First zero out S (has this already been done? don't think so)
  
  ZEROVEC(S, data->XS_blkptr[data->numblk+1]-1);

  // Do those A_i which are diagonal or sparse

  for(i = 1-obj; i <= data->m; i++) // Note 'i = 1-obj'
    for(j = data->AA_rowptr[i]; j <= data->AA_rowptr[i+1]-1; j++)
      S[data->AA_colind[j]] += y[i]*data->AA_colval_one[j];

  // Now take care of lowrank matrices (only when corresponding blocks of S are dense)

  for(h = 1; h <= data->lr_num; h++) {

    k  = data->lr_blk[h];  
    i  = data->lr_mat[h];  

    if(data->XS_blksto[k] == DENSE) {

      // sanity check (delete later)
      if(data->blktype[k] != SDPBLK) { printf("AToper: Unexpected block type!\n"); exit(0); }

      // Setup easy posize_ter to data matrix
      
      if(i == 0) A = data->C   [k];
      else       A = data->A[i][k];

      // Do the rank-1 updates of S
      
      for(j = 1; j <= A->lr->ncol; j++) {
        tv = y[i]*A->lr->d[j];
        dsyr_(&uplo, data->blksz+k, &tv, A->lr->ent + (j-1)*data->blksz[k] + 1, &inc1, S + data->XS_blkptr[k], data->blksz+k);
      }

    }

  }

  // If obj==1, make sure y[0] = origy0 
  
  if(obj) y[0] = origy0;

  return 1;
}

// size_t sparsesymmmat_timesmat(sparsesymmmat *S, double* vec, double* result, size_t n, size_t m)
// {
//   size_t k, r, c, nnz;
//   double ent;

//   // Note: this routine assumes first n elements of S are the diagonal

//   // Note: this routine adds onto the preexisting contents of "result"

//   nnz = S->nnz;

//   for(k = 1; k <= n; k++)
//     mydaxpy(m, S->ent[k], vec + S->row[k], n, result + S->col[k], n);


//   for(k = n+1; k <= nnz; k++) {

//     r = S->row[k]; c = S->col[k]; ent = S->ent[k];

//     mydaxpy(m, ent, vec + r, n, result + c, n);
//     mydaxpy(m, ent, vec + c, n, result + r, n);
//   }

//   return 1;

// }

double C_normdatamat(problemdata* data)
{ 
  size_t i, j, k;
  double Cnorm, *VD;
  lowrankmat *A;
  size_t temp, one=1;

  Cnorm = 0.0;
  for(k = 1; k <= data->numblk; k++) {
    if(data->blktype[k] == 's') {
      if(data->C[k]->type == 's') {
        temp = idamax_(&(data->C[k]->sp->nnz), data->C[k]->sp->ent + 1, &one);
        Cnorm = mymax( Cnorm, fabs(data->C[k]->sp->ent[temp]) );
      }
      else if(data->C[k]->type == 'l') {
        A = data->C[k]->lr;
        MYCALLOC(VD, double, A->nrow*A->ncol + 1);
        EASYDCOPY(A->nrow*A->ncol, A->ent, VD);
        for(j = 1; j <= A->ncol; j++)
          EASYDSCAL(A->nrow, A->d[j], VD + (j-1)*A->nrow);
        for(i = 1; i <= A->nrow; i++)
          for(j = i; j <= A->nrow; j++)
            Cnorm = mymax( Cnorm, fabs( myddot(A->ncol, A->ent + i, A->nrow, VD + j, A->nrow) ) );
        MYFREE(VD);
      }
    }
    else if(data->blktype[k] == 'd') {
      temp = idamax_(&(data->C[k]->diag->nnz), data->C[k]->diag->ent + 1, &one);
      Cnorm = mymax( Cnorm, fabs(data->C[k]->diag->ent[temp]) );
    }
  }

  return Cnorm;

}

double normdatamat(problemdata* data, size_t matnum)
{ 
  size_t i, j, k;
  char uplo='L', trans='T';
  double normsq, *VD, *VtV, *DVtVD, one=1.0, zero=0.0;
  datamat *DD;
  lowrankmat *A;
  sparsesymmmat *B;

  normsq = 0.0;

  for(k = 1; k <= data->numblk; k++) {

    if(matnum > 0) DD = data->A[matnum][k];
    else           DD = data->C        [k];

    if(data->blktype[k] == 's') {

      if(DD->type == 's') {

        B = DD->sp;

        for(j = 1; j <= B->nnz; j++)
          if(B->row[j] == B->col[j]) normsq +=     pow(B->ent[j], 2.0);
          else                       normsq += 2.0*pow(B->ent[j], 2.0);

      }
      else if(DD->type == 'l') {

        A = DD->lr;

        // Create VD (note C-style numbering)
        MYCALLOC(VD, double, A->nrow*A->ncol);
        EASYDCOPY(A->nrow*A->ncol, A->ent, VD-1);
        for(j = 0; j < A->ncol; j++)
          EASYDSCAL(A->nrow, A->d[j+1], VD + j*A->nrow - 1);
        // Create VtV (including upper triangle)
        MYCALLOC(VtV, double, A->ncol*A->ncol);
        dsyrk_(&uplo, &trans, &(A->ncol), &(A->nrow), &one, A->ent+1, &(A->nrow), &zero, VtV, &(A->ncol));
        for(i = 0; i < A->ncol; i++) for(j = i+1; j < A->ncol; j++) VtV[j*A->ncol + i] = VtV[i*A->ncol + j];
        // Create DVtVD (including upper triangle)
        MYCALLOC(DVtVD, double, A->ncol*A->ncol);
        dsyrk_(&uplo, &trans, &(A->ncol), &(A->nrow), &one, VD, &(A->nrow), &zero, DVtVD, &(A->ncol));
        for(i = 0; i < A->ncol; i++) for(j = i+1; j < A->ncol; j++) DVtVD[j*A->ncol + i] = DVtVD[i*A->ncol + j];
        // Dot them
        normsq += EASYDDOT(A->ncol*A->ncol, VtV-1, DVtVD-1);
        // Free memory
        MYFREE(VD); MYFREE(VtV); MYFREE(DVtVD);

      }

    }
    else if(data->blktype[k] == 'd') {
      for(j = 1; j <= DD->diag->nnz; j++)
        normsq += pow(DD->diag->ent[j], 2.0);
    }

  }

  return sqrt(normsq);

}

size_t essential_calcs(problemdata* data, double* R, double normC, double normb, double* val, double* rho_c_val, double* rho_f_val)
{
  *val = function(data, R);
  gradient(data, R);
  *rho_c_val = EASYDNRM2(data->nr, data->G)/(1.0 + normC);
  // *rho_c_val = EASYDNRM2(data->nr, data->G)/(1.0 + fabs(*val));
  *rho_f_val = EASYDNRM2(data->m, data->vio)/(1.0 + normb);

  return 0;
}

