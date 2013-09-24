#include "myinclude.h"

// #define mymax(A, B) ((A) > (B) ? (A) : (B))
// #define mymin(A, B) ((A) < (B) ? (A) : (B))


size_t copyscaledvectovec(double* dy, double da, double* dx, size_t n)
{
  EASYDCOPY(n, dx, dy);
  EASYDSCAL(n, da, dy);
  return 0;
}

size_t move_in_dir(double* dy, double* dx, double da, double* dir, size_t n)
{
  if(dy == dx) EASYDAXPY(n, da, dir, dy);
  else if(dy == dir) {
    EASYDSCAL(n, da, dy);
    EASYDAXPY(n, 1.0, dx, dy);
  }
  else {
    EASYDCOPY(n, dx, dy);
    EASYDAXPY(n, da, dir, dy);
  }
  return 1;
}

size_t mydaxpy(size_t n, double da, double* dx, size_t incx, double* dy, size_t incy)
{
  return daxpy_(&n,&da,dx,&incx,dy,&incy);
}

size_t mydcopy(size_t n, double* dx, size_t incx, double* dy, size_t incy)
{
  return dcopy_(&n,dx,&incx,dy,&incy);
}

double myddot(size_t n, double* dx, size_t incx, double* dy, size_t incy)
{
  return ddot_(&n, dx, &incx, dy, &incy);
}

double mydnrm2(size_t n, double* dx, size_t incx)
{
  return dnrm2_(&n, dx, &incx);
}

size_t mydscal(size_t n, double da, double* dx, size_t incx)
{
  return dscal_(&n, &da, dx, &incx);
}


size_t createlowrankmat(lowrankmat** passedR, size_t ncol, size_t nrow)
{
  lowrankmat *R;

  MYCALLOC(R, lowrankmat, 1);
  
  R->ncol = ncol;
  R->nrow = nrow;
  
  MYCALLOC(R->d, double, ncol + 1);
  MYCALLOC(R->ent, double, nrow*ncol + 1);

  *passedR = R;

  return 1;
}

size_t destroylowrankmat(lowrankmat* R)
{
  MYFREE(R->d);
  MYFREE(R->ent);
  MYFREE(R);

  return 1;
}

size_t createsparsesymmmat(sparsesymmmat** passedS, size_t nnz)
{
  sparsesymmmat *S;

  MYCALLOC(S, sparsesymmmat, 1);
  MYCALLOC(S->row, size_t, nnz+1);
  MYCALLOC(S->col, size_t, nnz+1);
  S->nnz = nnz;
  MYCALLOC(S->ent, double, nnz+1);
  MYCALLOC(S->XS_in, size_t, nnz+1);

  *passedS = S;

  return 1;

}

size_t destroysparsesymmmat(sparsesymmmat* S)
{
  MYFREE(S->row);
  MYFREE(S->col);
  MYFREE(S->ent);
  MYFREE(S->XS_in);
  MYFREE(S);

  return 1;
}

size_t creatediagmat(diagmat** passedD, size_t nnz)
{
  diagmat *D;

  MYCALLOC(D, diagmat, 1);
  MYCALLOC(D->ind, size_t, nnz+1);
  D->nnz = nnz;
  MYCALLOC(D->ent, double, nnz+1);
  MYCALLOC(D->XS_in, size_t, nnz+1);

  *passedD = D;

  return 1;

}

size_t destroydiagmat(diagmat* D)
{
  MYFREE(D->ind);
  MYFREE(D->ent);
  MYFREE(D->XS_in);
  MYFREE(D);

  return 1;
}


size_t createdatamat(datamat** passedA, char type, size_t ncol_or_nnz, size_t dim, char* label)
{
  datamat *A;

  MYCALLOC(A, datamat, 1);
  A->type = type;
  MYCALLOC(A->label, char, 30);
  strcpy(A->label, label);

  if(type == 'l')
    createlowrankmat(&(A->lr), ncol_or_nnz, dim);
  
  if(type == 's')
    createsparsesymmmat(&(A->sp), ncol_or_nnz);

  if(type == 'd')
    creatediagmat(&(A->diag), ncol_or_nnz);

  // if type = 'u' then do nothing

  *passedA = A;

  return 1;
}

size_t destroydatamat(datamat* A)
{
  if(A->type == 'l')
    destroylowrankmat(A->lr);
  
  if(A->type == 's')
    destroysparsesymmmat(A->sp);

  if(A->type == 'd')
    destroydiagmat(A->diag);

  MYFREE(A->label);
  MYFREE(A);

  return 1;
}

