#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "mex.h"
#include "matrix.h"

#define MATIND(I,J,DIM) (J-1)*DIM + I - 1

#ifdef __WIN32
#define dsyrk_ dsyrk
#endif

#define MYCALLOC(VAR,TYPE,SIZE) VAR = (TYPE*)calloc(SIZE, sizeof(TYPE))
#define MYFREE(VAR) free(VAR)

#define mymin(A, B) ((A) < (B) ? (A) : (B))
#define mymax(A, B) ((A) > (B) ? (A) : (B))

#define DATABLOCKIND(DATA,BLOCK,NUMBLOCK) ((DATA+1)-1)*NUMBLOCK + BLOCK - 1

extern void dsyrk_(); 
size_t classifyApq(size_t*, char*, size_t*, size_t*, size_t, size_t, size_t, size_t*, size_t*, size_t*, size_t*);
size_t getparams(char *paramfile, size_t *inputtype, double *rho_f, double *rho_c, double *sigmafac, size_t *rankreduce, size_t *timelim, size_t *printlevel, size_t *dthresh_dim, double *dthresh_dens, size_t *numbfgsvecs, double *rankredtol, double *gaptol, ptrdiff_t *checkbd, size_t *typebd);
size_t getstorage(size_t m, size_t numblk, size_t* blksz, char* blktype, 
               size_t* CAinfo_entptr, size_t* passedn, size_t* passednr,
               size_t* maxranks);
size_t sdplrlib(size_t m, size_t numblk, size_t *blksz, char *blktype, double *b, double *CAent, size_t *CArow, size_t *CAcol, size_t *CAinfo_entptr, char *CAinfo_type, size_t numbfgsvecs, double rho_f, double rho_c, double sigmafac, size_t rankreduce, double gaptol, ptrdiff_t checkbd, size_t typebd, size_t dthresh_dim, double dthresh_dens, size_t timelim, double rankredtol, size_t printlevel, double *R, double *lambda, size_t *maxranks, size_t *ranks, double *pieces);
size_t writedata_raw(char*, size_t, size_t, size_t*, char*, double*, double*, size_t*, size_t*, size_t*, size_t*, char*, char*);

size_t writedata_sdpa(char* datafilename, size_t m, size_t numblk, size_t* blksz,
                   char* blktype, double* b, double* CAent,
                   size_t* CArow, size_t* CAcol, size_t* CAinfo_entptr,
                   char* CAinfo_type);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  /* Declarations */

  size_t m;
  size_t numblk;
  size_t *blksz;
  char *blktype;
  double *b;
  double *CAent;
  size_t *CArow;
  size_t *CAcol;
  size_t *CAinfo_entptr;
  size_t *CAinfo_rowcolptr;
  char *CAinfo_type;
  char *CAinfo_storage;

  size_t *maxranks;
  size_t *ranks;
  size_t nr;

  double *R;
  double *lambda;
  double pieces[8];

  size_t numbfgsvecs;
  double rho_f;
  double rho_c;
  double sigmafac;
  size_t reorder;
  size_t precond;
  double gdens;
  size_t rankreduce;
  double gaptol;
  ptrdiff_t checkbd;
  size_t typebd;
  size_t dthresh_dim;
  double dthresh_dens;
  size_t timelim;
  double rankredtol;
  size_t inputtype;
  size_t doAR;
  size_t printlevel;

  size_t writeraw=0;
  size_t soln_factored=0;
  double *forcerank;
  size_t seed;

  size_t h;
  size_t i;
  size_t j;
  size_t k;

  size_t n;

  const mxArray *A_;
  const mxArray *b_;
  const mxArray *c_;
  const mxArray *K_;
  const mxArray *pars_;

  const mxArray *lrA_ = NULL;

  const mxArray *x0_;
  const mxArray *y0_;
  const mxArray *info0_;
  const mxArray *r0_;

  const mxArray *Kl_;
  const mxArray *Ks_;
  size_t Ks_first=-1;
  size_t n_;


  size_t A_sp;
  size_t A_tr;
  size_t b_sp;
  size_t b_tr;
  size_t c_sp;
  size_t c_tr;

  size_t c_empty;
  size_t pars_empty=1;
  size_t lrA_empty=1;
  size_t y0_empty=1;
  size_t info0_empty=1;
  size_t r0_empty=1;

  mxArray *field;
  mxArray *cell;
  char *fieldname;
  double *mat;
  size_t nnz;
  size_t *ir;
  size_t *jc;
  size_t *inblk;
  size_t *inblkptr;
  size_t cons;
  size_t blk;
  size_t ii;
  size_t jj;
  size_t *tempsize_t;
  size_t base;
  size_t base1;
  char uplo = 'l';
  char trans = 'n';
  double one = 1.0;
  double zero = 0.0;

  char **infonames;

  /*
     Basic error checking
  */

  /* correct # of outputs */
  if(nlhs > 4)
    mexErrMsgTxt("sdplr: At most four output arguments.");

  /* correct # of inputs */
  if(nrhs < 4 || nrhs > 10)
    mexErrMsgTxt("sdplr: Four to ten arguments.");

  A_ = prhs[0];
  b_ = prhs[1];
  c_ = prhs[2];

  m = mymin( mxGetM(A_) , mxGetN(A_) );
  n_ = mymax( mxGetM(A_) , mxGetN(A_) );

  /* matrix passed in for A (not complex) */
  if( (mxIsSparse(A_) == 1 && mxGetPi(A_) != NULL) ||
      (mxIsSparse(A_) == 0 && (mxIsDouble(A_) == 0 || mxIsComplex(A_) == 1)) )
    mexErrMsgTxt("sdplr: Input 1 must be of type double.");

  /* matrix passed in for b (not complex) */
  if( (mxIsSparse(b_) == 1 && mxGetPi(b_) != NULL) ||
      (mxIsSparse(b_) == 0 && (mxIsDouble(b_) == 0 || mxIsComplex(b_) == 1)) )
    mexErrMsgTxt("sdplr: Input 2 must be of type double .");

  /* b is a col vec */
  if( mymin( mxGetM(b_) , mxGetN(b_) ) != 1)
    mexErrMsgTxt("sdplr: Input 2 must be a vector.");

  /* small dim of A equals big dim of b */
  if( m != mymax( mxGetM(b_) , mxGetN(b_) ) )
    mexErrMsgTxt("sdplr: Inputs 1 and 2 do not have compatible sizes.");

  /* matrix passed in for c (not complex) */
  if( (mxIsSparse(c_) == 1 && mxGetPi(c_) != NULL) ||
      (mxIsSparse(c_) == 0 && (mxIsDouble(c_) == 0 || mxIsComplex(c_) == 1)) )
    mexErrMsgTxt("sdplr: Input 3 must be of type double.");

  /* c is a col vec */
  if( mymin( mxGetM(c_) , mxGetN(c_) ) == 0 )
    c_empty = 1;
  else {
    if( mymin( mxGetM(c_) , mxGetN(c_) ) != 1 )
      mexErrMsgTxt("sdplr: Input 3 must be a vector.");
    c_empty = 0;
  }

  /* big dim of A equals big dim of c */
  if( c_empty == 0 && n_ != mymax( mxGetM(c_) , mxGetN(c_) ) )
    mexErrMsgTxt("sdplr: Inputs 1 and 3 do not have compatible sizes.");

  if(nrhs > 3) {

    K_ = prhs[3];
    Kl_ = NULL;
    Ks_ = NULL;

    /* structure passed for K */
    if(mxIsStruct(K_) == 0)
      mexErrMsgTxt("sdplr: Input 4 must be of type struct.");

    /* K must have only one set of information */
    if(mxGetM(K_)*mxGetN(K_) != 1)
      mexErrMsgTxt("sdplr: Input 4 must be a structure of size 1x1.");

    /* Only linear and semidefinite cones specified */
    for(i = 0; i < mxGetNumberOfFields(K_); i++) {
      field = mxGetFieldByNumber(K_, 0, i);
      if( mxGetM(field)*mxGetN(field) != 0 ) {

        mat = mxGetPr(field);
        if(mat[0] != 0.0) {

          fieldname = (char*)mxGetFieldNameByNumber(K_, i); /* Don't know why mxGetFieldNameByNumber is "const" */
          if(strcmp(fieldname, "s") != 0 && strcmp(fieldname, "l") != 0)
            mexErrMsgTxt("sdplr: Only linear and semidefinite cones supported at this time.");
          if(strcmp(fieldname, "l") == 0) {
            Kl_ = field;
            h = i;
            if(Ks_first == -1) Ks_first = 0;
          }
          if(strcmp(fieldname, "s") == 0) {
            Ks_ = field;
            j = i;
            if(Ks_first == -1) Ks_first = 1;
          }

        }

      }
    }

    /* Linear must be specified before psd
    if(Kl_ != NULL && Ks_ != NULL && h > j)
      mexErrMsgTxt("sdplr: Linear portion of cone must be specified before semidefinite portion."); */

    /* Portions of cone must be double matrices (not complex) */
    if(Kl_ != NULL)
      if(mxIsDouble(Kl_) == 0 || mxIsComplex(Kl_) == 1)
        mexErrMsgTxt("sdplr: Linear portion of cone must be of type double.");
    if(Ks_ != NULL)
      if(mxIsDouble(Ks_) == 0 || mxIsComplex(Ks_) == 1)
        mexErrMsgTxt("sdplr: Semidefinite portion of cone must be of type double.");

    /* Linear portion of K must be scalar and psd portion must be vector */
    if(Kl_ != NULL)
      if( mxGetM(Kl_)*mxGetN(Kl_) != 1 )
        mexErrMsgTxt("sdplr: Linear portion of cone must be scalar indicating size.");
    if(Ks_ != NULL)
      if( mymin( mxGetM(Ks_) , mxGetN(Ks_) ) != 1 )
        mexErrMsgTxt("sdplr: Semidefinite portion of cone must be vector of sizes.");

    /* Total dim of K must equal max dim of A */
    h = 0;
    if(Kl_ != NULL) {
      mat = mxGetPr(Kl_);
      h += (size_t)mat[0];
    }
    if(Ks_ != NULL) {
      mat = mxGetPr(Ks_);
      for(i = 0; i < mymax( mxGetM(Ks_) , mxGetN(Ks_) ); i++)
        h += (size_t)mat[i]*(size_t)mat[i];
    }
    if(h != n_)
      mexErrMsgTxt("sdplr: Inputs 1 and 4 have incompatible sizes.");


    if(nrhs > 4) {

      pars_ = prhs[4];

      /* first determine if pars is empty */
      if( mymin( mxGetM(pars_) , mxGetN(pars_) ) == 0 )
        pars_empty = 1;
      else
        pars_empty = 0;

      if(!pars_empty) {

        /* structure passed for pars */
        if(mxIsStruct(pars_) == 0)
          mexErrMsgTxt("sdplr: Input 5 must be of type structure.");

        if( mxGetM(pars_) * mxGetN(pars_) != 1 )
          mexErrMsgTxt("sdplr: Input 5 must contain only one set of information.");

      }

      if(nrhs > 5) {

        lrA_ = prhs[5];
            
        /* first determine if lrA is empty */
        if( mymin( mxGetM(lrA_) , mxGetN(lrA_) ) == 0 )
          lrA_empty = 1;
        else
          lrA_empty = 0;

        if(!lrA_empty) {

          if(mxIsStruct(lrA_) == 0)
            mexErrMsgTxt("sdplr: Input 6 must be of type struct.");

          for(i = 0; i < mxGetM(lrA_)*mxGetN(lrA_); i++) {

            field = mxGetField(lrA_,i,"cons");
            if(field == NULL || mxGetM(field)*mxGetN(field) != 1 ||
              mxIsDouble(field) == 0 || mxIsComplex(field) == 1 ||
              (size_t)mxGetScalar(field) < 0 || (size_t)mxGetScalar(field) > m)
              mexErrMsgTxt("sdplr: Input 6 has incorrect format.");

            field = mxGetField(lrA_,i,"start");
            if(field == NULL || mxGetM(field)*mxGetN(field) != 1 ||
              mxIsDouble(field) == 0 || mxIsComplex(field) == 1 ||
              (size_t)mxGetScalar(field) < 0 || (size_t)mxGetScalar(field) > n_)
              mexErrMsgTxt("sdplr: Input 6 has incorrect format.");

            field = mxGetField(lrA_,i,"D");
            if(field == NULL || mymin( mxGetM(field) , mxGetN(field) ) != 1 ||
              mxIsDouble(field) == 0 || mxIsComplex(field) == 1)
              mexErrMsgTxt("sdplr: Input 6 has incorrect format.");

            h = mymax( mxGetM(field) , mxGetN(field) );

            field = mxGetField(lrA_,i,"V");
            if(field == NULL || mxGetN(field) != h ||
              mxIsDouble(field) == 0 || mxIsComplex(field) == 1)
              mexErrMsgTxt("sdplr: Input 6 has incorrect format.");

          }        

        }

        if(nrhs > 6) {

          /* User has input x0 */
          /* Because of pre-processing in sdplr.m, x0 is 0x0 */
          /* Double check this */
          
          x0_ = prhs[6];
          if( mxGetM(x0_) > 0 || mxGetN(x0_) > 0 )
            mexErrMsgTxt("sdplr (size_ternal error): Input 7 has not been replaced correctly.");

          if(nrhs > 7) {

            /* User has input y0 */
            
            y0_ = prhs[7];

            /* first determine if y0 is empty */
            if( mymin( mxGetM(y0_) , mxGetN(y0_) ) == 0 )
              y0_empty = 1;
            else
              y0_empty = 0;

            if(!y0_empty) {

              /* matrix passed in for y0 (not complex) */
              if( (mxIsSparse(y0_) == 1 && mxGetPi(y0_) != NULL) ||
                  (mxIsSparse(y0_) == 0 && (mxIsDouble(y0_) == 0 || mxIsComplex(y0_) == 1)) )
                mexErrMsgTxt("sdplr: Input 8 must be of type double.");

              /* y0 is column vector of size m */
              if( mymin( mxGetM(y0_) , mxGetN(y0_) ) != 1 && mymax( mxGetM(y0_) , mxGetN(y0_) ) != m )
                mexErrMsgTxt("sdplr: Input 8 must be a vector with the same length as input 2.");

            }

            if(nrhs > 8) {

              /* User has input info0 */
            
              info0_ = prhs[8];

              /* determine if info0 is empty */
              if( mymin( mxGetM(info0_) , mxGetN(info0_) ) == 0 )
                info0_empty = 1;
              else
                info0_empty = 0;

              if(nrhs > 9) {

                /* User has input r0 */
            
                r0_ = prhs[9];

                /* determine if r0 is empty */
                if( mymin( mxGetM(r0_) , mxGetN(r0_) ) == 0 )
                  r0_empty = 1;
                else
                  r0_empty = 0;

              }

            }

          }

        }

      }

    }

  }

  /* Create some helpful flags */
  if(mxIsSparse(A_) == 1) A_sp = 1;
  else A_sp = 0;
  if( mxGetM(A_) == m ) A_tr = 0;
  else A_tr = 1;
  if(mxIsSparse(b_) == 1) b_sp = 1;
  else b_sp = 0;
  if( mxGetM(b_) == m ) b_tr = 0;
  else b_tr = 1;
  if(c_empty == 0) {
    if(mxIsSparse(c_) == 1) c_sp = 1;
    else c_sp = 0;
    if(mxGetM(c_) == n_ ) c_tr = 0;
    else c_tr = 1;
  }


  /*
   Create problem structures (phase 1)
   */

  /* Already have gotten m */

  /* Get numblk */
  numblk = 0;
  if(Kl_ != NULL) numblk++;
  if(Ks_ != NULL) numblk += mymax( mxGetM(Ks_) , mxGetN(Ks_) );

  /*
   Setup structure to classify entries of A
   */
  MYCALLOC(inblk, size_t, n_);
  MYCALLOC(inblkptr, size_t, numblk+2);
  h = 0; k = 1;
  /* This code is repetitive; just trying to fix a bug regarding order of cones */
  if(Ks_first == 1) {
    if(Ks_ != NULL) {
      mat = mxGetPr(Ks_);
      for(i = 0; i < mymax( mxGetM(Ks_) , mxGetN(Ks_) ); i++) {
        inblkptr[k] = h;
        for(j = 0; j < (size_t)mat[i]*(size_t)mat[i]; j++) {
          inblk[h] = k;
          h++;
        }
        k++;
      }
    }
    if(Kl_ != NULL) {
      mat = mxGetPr(Kl_);
      inblkptr[k] = h;
      for(j = 0; j < (size_t)mat[0]; j++) {
        inblk[h] = k;
        h++;
      }
      k++;
    }
  }
  else {
    if(Kl_ != NULL) {
      mat = mxGetPr(Kl_);
      inblkptr[k] = h;
      for(j = 0; j < (size_t)mat[0]; j++) {
        inblk[h] = k;
        h++;
      }
      k++;
    }
    if(Ks_ != NULL) {
      mat = mxGetPr(Ks_);
      for(i = 0; i < mymax( mxGetM(Ks_) , mxGetN(Ks_) ); i++) {
        inblkptr[k] = h;
        for(j = 0; j < (size_t)mat[i]*(size_t)mat[i]; j++) {
          inblk[h] = k;
          h++;
        }
        k++;
      }
    }
  }
  inblkptr[k] = h;
  if(h != n_ || k-1 != numblk)
    mexErrMsgTxt("sdplr: Internal error (4).");

  /*
   Create problem structures (phase 2)
   */

  /* Get blksz and blktype */
  MYCALLOC(blksz, size_t, numblk);
  MYCALLOC(blktype, char, numblk);
  h = 0;
  /* This code is repetitive; just trying to fix a bug regarding order of cones */
  if(Ks_first == 1) {
    if(Ks_ != NULL) {
      mat = mxGetPr(Ks_);
      for(i = 0; i < mymax( mxGetM(Ks_) , mxGetN(Ks_) ); i++) {
        blksz[h] = (size_t)mat[i];
        blktype[h] = 's';
        h++;
      }
    }
    if(Kl_ != NULL) {
      mat = mxGetPr(Kl_);
      blksz[h] = (size_t)mat[0];
      blktype[h] = 'd';
      h++;
    }
  }
  else {
    if(Kl_ != NULL) {
      mat = mxGetPr(Kl_);
      blksz[h] = (size_t)mat[0];
      blktype[h] = 'd';
      h++;
    }
    if(Ks_ != NULL) {
      mat = mxGetPr(Ks_);
      for(i = 0; i < mymax( mxGetM(Ks_) , mxGetN(Ks_) ); i++) {
        blksz[h] = (size_t)mat[i];
        blktype[h] = 's';
        h++;
      }
    }
  }
  if(h != numblk)
    mexErrMsgTxt("sdplr: Internal error (3).");

  /* Get b */
  MYCALLOC(b, double, m);
  if(b_sp) {
    mat = mxGetPr(b_);
    if(b_tr) {
      jc = mxGetJc(b_);
      for(j = 0; j < m; j++)
        for(i = jc[j]; i <= jc[j+1]-1; i++)
          b[j] = mat[i];
    }
    else {
      nnz = mxGetNzmax(b_);
      ir = mxGetIr(b_);
      for(i = 0; i < nnz; i++)
        b[ ir[i] ] = mat[i];
    }
  }
  else {
    mat = mxGetPr(b_);
    for(i = 0; i < m; i++)
      b[i] = mat[i];
  }

  /* Count nnz in CAent, etc. */
  nnz = 0;
  if(c_empty == 0) {
    if(c_sp) {
      jc = mxGetJc(c_);
      ir = mxGetIr(c_);
      for(j = 0; j < mxGetN(c_); j++)
        for(i = jc[j]; i <= jc[j+1]-1; i++) {
          classifyApq(blksz, blktype, inblk, inblkptr, ir[i], j, 1 - c_tr, &cons, &blk, &ii, &jj);
          if(ii >= jj) nnz++;
        }
    }
    else {
      for(i = 0; i < n_; i++) {
        classifyApq(blksz, blktype, inblk, inblkptr, 0, i, 0, &cons, &blk, &ii, &jj);
        if(ii >= jj) nnz++;
      }
    }
  }
  if(A_sp) {
    jc = mxGetJc(A_);
    ir = mxGetIr(A_);
    for(j = 0; j < mxGetN(A_); j++)
      for(i = jc[j]; (int)i <= (int)jc[j+1]-(int)1; i++) {
        classifyApq(blksz, blktype, inblk, inblkptr, ir[i], j, A_tr, &cons, &blk, &ii, &jj);
        if(ii >= jj) nnz++;
      }
  }
  else {
    for(j = 0; j < mxGetN(A_); j++)
      for(i = 0; i < mxGetM(A_); i++) {
        classifyApq(blksz, blktype, inblk, inblkptr, i, j, A_tr, &cons, &blk, &ii, &jj);
        if(ii >= jj) nnz++;
      }
  }
  if(lrA_ != NULL && !lrA_empty) {

    for(h = 0; h < mxGetM(lrA_)*mxGetN(lrA_); h++) {

      field = mxGetField(lrA_,h,"start");
      if(A_tr) i = (size_t)mxGetScalar(field) - 1;
      else     j = (size_t)mxGetScalar(field) - 1;

      field = mxGetField(lrA_,h,"cons");
      if(A_tr) j = (size_t)mxGetScalar(field) - 1;
      else     i = (size_t)mxGetScalar(field) - 1;

      classifyApq(blksz, blktype, inblk, inblkptr, i, j, A_tr, &cons, &blk, &ii, &jj);

      if(ii != 1 || jj != 1)
        mexErrMsgTxt("sdplr: Internal error (6).");
      if(cons != (size_t)mxGetScalar(field))
        mexErrMsgTxt("sdplr: Internal error (6).");

      field = mxGetField(lrA_,h,"D");
      k = mymax( mxGetM(field) , mxGetN(field) );

      field = mxGetField(lrA_,h,"V");
      if( blktype[blk-1] != 's' || mxGetM(field) != blksz[blk-1] )
        mexErrMsgTxt("sdplr: Internal error (7).");

      nnz += k*(blksz[blk-1] + 1);
      
    }

  }


  /* Allocate CA space */
  MYCALLOC(CAent, double, nnz);
  MYCALLOC(CArow, size_t, nnz);
  MYCALLOC(CAcol, size_t, nnz);

  MYCALLOC(CAinfo_entptr, size_t, (m+1)*numblk + 1 );
  MYCALLOC(CAinfo_rowcolptr, size_t, (m+1)*numblk + 1 );
  MYCALLOC(CAinfo_type, char, (m+1)*numblk );
  MYCALLOC(CAinfo_storage, char, (m+1)*numblk );

  if( DATABLOCKIND(m,numblk,numblk) != (m+1)*numblk - 1 )
    mexErrMsgTxt("sdplr: Internal error (2).");

  /* Easy: setup CAinfo_storage (default) */
  for(i = 0; i < (m+1)*numblk; i++) CAinfo_storage[i] = 's';

  /* Easy: setup CAinfo_type (default) */
  for(i = 0; i <= m; i++)
    for(k = 1; k <= numblk; k++) {
      if(blktype[k-1] == 'd')
        CAinfo_type[ DATABLOCKIND(i,k,numblk) ] = 'd';
      else
        CAinfo_type[ DATABLOCKIND(i,k,numblk) ] = 's'; /* will have low-rank here in the future */
    }

  /* Adjust CAinfo_storage and CAinfo_type for low-rank matrices (kind of easy) */
  if(lrA_ != NULL && !lrA_empty) {

    for(h = 0; h < mxGetM(lrA_)*mxGetN(lrA_); h++) {

      field = mxGetField(lrA_,h,"start");
      if(A_tr) i = (size_t)mxGetScalar(field) - 1;
      else     j = (size_t)mxGetScalar(field) - 1;

      field = mxGetField(lrA_,h,"cons");
      if(A_tr) j = (size_t)mxGetScalar(field) - 1;
      else     i = (size_t)mxGetScalar(field) - 1;

      classifyApq(blksz, blktype, inblk, inblkptr, i, j, A_tr, &cons, &blk, &ii, &jj);

      CAinfo_storage[ DATABLOCKIND(cons,blk,numblk) ] = 'd';
      CAinfo_type[ DATABLOCKIND(cons,blk,numblk) ] = 'l';
    }

  }

  /* Setup CAinfo_entptr and CAinfo_rowcolptr */

  MYCALLOC(tempsize_t, size_t, (m+1)*numblk);

  if(c_empty == 0) {
    if(c_sp) {
      jc = mxGetJc(c_);
      ir = mxGetIr(c_);
      for(j = 0; j < mxGetN(c_); j++)
        for(i = jc[j]; i <= jc[j+1]-1; i++) {
          classifyApq(blksz, blktype, inblk, inblkptr, ir[i], j, 1 - c_tr, &cons, &blk, &ii, &jj);
          cons = 0;
          if(ii >= jj) tempsize_t[ DATABLOCKIND(cons,blk,numblk) ]++;
        }
    }
    else{
      for(i = 0; i < n_; i++) {
        classifyApq(blksz, blktype, inblk, inblkptr, 0, i, 0, &cons, &blk, &ii, &jj);
        cons = 0;
        if(ii >= jj) tempsize_t[ DATABLOCKIND(cons,blk,numblk) ]++;
      }
    }
  }
  if(A_sp) {
    jc = mxGetJc(A_);
    ir = mxGetIr(A_);
    for(j = 0; j < mxGetN(A_); j++)
      for(i = jc[j]; (int)i <= (int)jc[j+1]-(int)1; i++) {
        classifyApq(blksz, blktype, inblk, inblkptr, ir[i], j, A_tr, &cons, &blk, &ii, &jj);
        if(ii >= jj) tempsize_t[ DATABLOCKIND(cons,blk,numblk) ]++;
      }
  }
  else {
    for(j = 0; j < mxGetN(A_); j++)
      for(i = 0; i < mxGetM(A_); i++) {
        classifyApq(blksz, blktype, inblk, inblkptr, i, j, A_tr, &cons, &blk, &ii, &jj);
        if(ii >= jj) tempsize_t[ DATABLOCKIND(cons,blk,numblk) ]++;
      }
  }
  if(lrA_ != NULL && !lrA_empty) {

    for(h = 0; h < mxGetM(lrA_)*mxGetN(lrA_); h++) {

      field = mxGetField(lrA_,h,"start");
      if(A_tr) i = (size_t)mxGetScalar(field) - 1;
      else     j = (size_t)mxGetScalar(field) - 1;

      field = mxGetField(lrA_,h,"cons");
      if(A_tr) j = (size_t)mxGetScalar(field) - 1;
      else     i = (size_t)mxGetScalar(field) - 1;

      classifyApq(blksz, blktype, inblk, inblkptr, i, j, A_tr, &cons, &blk, &ii, &jj);

      field = mxGetField(lrA_,h,"D");
      k = mymax( mxGetM(field) , mxGetN(field) );

      tempsize_t[ DATABLOCKIND(cons,blk,numblk) ] += k*(blksz[blk-1] + 1);
      
    }

  }


  h = 0; /* begin sanity check */
  for(i = 0; i < (m+1)*numblk; i++)
    h += tempsize_t[i];
  if(h != nnz)
    mexErrMsgTxt("sdplr: Internal error (5)."); /* end sanity check */

  h = 0;
  for(i = 0; i <= (m+1)*numblk; i++) {
    CAinfo_entptr[i] = h;
    CAinfo_rowcolptr[i] = h;
    if(i < (m+1)*numblk) h += tempsize_t[i];
  }
//   printf("CAinfo_entptr[0] = %d\n", CAinfo_entptr[0]);

  MYFREE(tempsize_t);

  /* Setup CAent, CArow, and CAcol */

  MYCALLOC(tempsize_t, size_t, (m+1)*numblk);

  if(c_empty == 0) {
    if(c_sp) {
      mat = mxGetPr(c_);
      jc = mxGetJc(c_);
      ir = mxGetIr(c_);
      for(j = 0; j < mxGetN(c_); j++)
        for(i = jc[j]; i <= jc[j+1]-1; i++) {
          classifyApq(blksz, blktype, inblk, inblkptr, ir[i], j, 1 - c_tr, &cons, &blk, &ii, &jj);
          cons = 0;
          if(ii >= jj) {
            h = CAinfo_entptr[ DATABLOCKIND(cons,blk,numblk) ] + tempsize_t[ DATABLOCKIND(cons,blk,numblk) ];
            CAent[h] = mat[i];
            h = CAinfo_rowcolptr[ DATABLOCKIND(cons,blk,numblk) ] + tempsize_t[ DATABLOCKIND(cons,blk,numblk) ];
            CArow[h] = ii;
            CAcol[h] = jj;
            tempsize_t[ DATABLOCKIND(cons,blk,numblk) ]++;
          }
        }
    }
    else{
      mat = mxGetPr(c_);
      for(i = 0; i < n_; i++) {
        classifyApq(blksz, blktype, inblk, inblkptr, 0, i, 0, &cons, &blk, &ii, &jj);
        cons = 0;
        if(ii >= jj) {
          h = CAinfo_entptr[ DATABLOCKIND(cons,blk,numblk) ] + tempsize_t[ DATABLOCKIND(cons,blk,numblk) ];
          CAent[h] = mat[i];
          h = CAinfo_rowcolptr[ DATABLOCKIND(cons,blk,numblk) ] + tempsize_t[ DATABLOCKIND(cons,blk,numblk) ];
          CArow[h] = ii;
          CAcol[h] = jj;
          tempsize_t[ DATABLOCKIND(cons,blk,numblk) ]++;
        }
      }
    }
  }
  if(A_sp) {
    mat = mxGetPr(A_);
    jc = mxGetJc(A_);
    ir = mxGetIr(A_);
    for(j = 0; j < mxGetN(A_); j++)
      for(i = jc[j]; (int)i <= (int)jc[j+1]-(int)1; i++) {
        classifyApq(blksz, blktype, inblk, inblkptr, ir[i], j, A_tr, &cons, &blk, &ii, &jj);
        if(ii >= jj) {
          h = CAinfo_entptr[ DATABLOCKIND(cons,blk,numblk) ] + tempsize_t[ DATABLOCKIND(cons,blk,numblk) ];
          CAent[h] = mat[i];
/*          printf("CAent[%d] = %f\n", h, CAent[h]); */
          h = CAinfo_rowcolptr[ DATABLOCKIND(cons,blk,numblk) ] + tempsize_t[ DATABLOCKIND(cons,blk,numblk) ];
          CArow[h] = ii;
          CAcol[h] = jj;
          tempsize_t[ DATABLOCKIND(cons,blk,numblk) ]++;
        }
      }
  }
  else {
    mat = mxGetPr(A_);
    for(j = 0; j < mxGetN(A_); j++)
      for(i = 0; i < mxGetM(A_); i++) {
        classifyApq(blksz, blktype, inblk, inblkptr, i, j, A_tr, &cons, &blk, &ii, &jj);
        if(ii >= jj) {
          h = CAinfo_entptr[ DATABLOCKIND(cons,blk,numblk) ] + tempsize_t[ DATABLOCKIND(cons,blk,numblk) ];
          CAent[h] = mat[j*mxGetM(A_) + i];
          h = CAinfo_rowcolptr[ DATABLOCKIND(cons,blk,numblk) ] + tempsize_t[ DATABLOCKIND(cons,blk,numblk) ];
          CArow[h] = ii;
          CAcol[h] = jj;
          tempsize_t[ DATABLOCKIND(cons,blk,numblk) ]++;
        }
      }
  }
  if(lrA_ != NULL && !lrA_empty) {

    for(h = 0; h < mxGetM(lrA_)*mxGetN(lrA_); h++) {

      field = mxGetField(lrA_,h,"start");
      if(A_tr) i = (size_t)mxGetScalar(field) - 1;
      else     j = (size_t)mxGetScalar(field) - 1;

      field = mxGetField(lrA_,h,"cons");
      if(A_tr) j = (size_t)mxGetScalar(field) - 1;
      else     i = (size_t)mxGetScalar(field) - 1;

      classifyApq(blksz, blktype, inblk, inblkptr, i, j, A_tr, &cons, &blk, &ii, &jj);

      field = mxGetField(lrA_,h,"D");
      k = mymax( mxGetM(field) , mxGetN(field) );

      /* i and j are available */
      mat = mxGetPr(field);
      for(i = 0; i < k; i++)
        CAent[ CAinfo_entptr[ DATABLOCKIND(cons,blk,numblk) ] + i ] = mat[i];

      field = mxGetField(lrA_,h,"V");
      mat = mxGetPr(field);
      for(i = 0; i < k*blksz[blk-1]; i++)
        CAent[ CAinfo_entptr[ DATABLOCKIND(cons,blk,numblk) ] + k + i ] = mat[i];

      /* don't have to do anything for CArow and CAcol b/c stored as dense */

      tempsize_t[ DATABLOCKIND(cons,blk,numblk) ] += k*(blksz[blk-1] + 1);
      
    }

  }

  MYFREE(tempsize_t);

  /* Get default parameters */

  getparams(NULL, &inputtype, &rho_f, &rho_c,
            &sigmafac, &rankreduce, &timelim, &printlevel, &dthresh_dim,
            &dthresh_dens, &numbfgsvecs, &rankredtol, &gaptol, &checkbd, &typebd);

  /* including MEX specific parameters */
  writeraw = 0;
  soln_factored = 0;
  forcerank = NULL;
  seed = (unsigned)time(NULL);


  /* Overwrite with params specified by user */

/*
feastol
centol
dir
penfac
reduce
limit
printlevel
dthresh_dim
dthresh_den
numlbfgs
feasswitch
auxtn
ranktol
precond
fillden
reorder
gaptol
checkbd
bdtype
seed
*/

  if(nrhs > 4 && !pars_empty) {

    /* pars_ already allocated */

    field = mxGetField(pars_, 0, "feastol");
    if( field != NULL && mxIsDouble(field) == 1 && mxIsComplex(field) == 0 && mxGetM(field)*mxGetN(field) == 1 )
      rho_f = mxGetScalar(field);

    field = mxGetField(pars_, 0, "centol");
    if( field != NULL && mxIsDouble(field) == 1 && mxIsComplex(field) == 0 && mxGetM(field)*mxGetN(field) == 1 )
      rho_c = mxGetScalar(field);

/*     field = mxGetField(pars_, 0, "dir");
     if( field != NULL && mxIsDouble(field) == 1 && mxIsComplex(field) == 0 && mxGetM(field)*mxGetN(field) == 1 )
       method = (size_t)mxGetScalar(field);
       */

    field = mxGetField(pars_, 0, "penfac");
    if( field != NULL && mxIsDouble(field) == 1 && mxIsComplex(field) == 0 && mxGetM(field)*mxGetN(field) == 1 )
      sigmafac = mxGetScalar(field);

    field = mxGetField(pars_, 0, "reduce");
    if( field != NULL && mxIsDouble(field) == 1 && mxIsComplex(field) == 0 && mxGetM(field)*mxGetN(field) == 1 )
      rankreduce = (size_t)mxGetScalar(field);

    field = mxGetField(pars_, 0, "limit");
    if( field != NULL && mxIsDouble(field) == 1 && mxIsComplex(field) == 0 && mxGetM(field)*mxGetN(field) == 1 )
      timelim = (size_t)mxGetScalar(field);

    field = mxGetField(pars_, 0, "printlevel");
    if( field != NULL && mxIsDouble(field) == 1 && mxIsComplex(field) == 0 && mxGetM(field)*mxGetN(field) == 1 )
      printlevel = (size_t)mxGetScalar(field);

    field = mxGetField(pars_, 0, "dthresh_dim");
    if( field != NULL && mxIsDouble(field) == 1 && mxIsComplex(field) == 0 && mxGetM(field)*mxGetN(field) == 1 )
      dthresh_dim = (size_t)mxGetScalar(field);

    field = mxGetField(pars_, 0, "dthresh_den");
    if( field != NULL && mxIsDouble(field) == 1 && mxIsComplex(field) == 0 && mxGetM(field)*mxGetN(field) == 1 )
      dthresh_dens = mxGetScalar(field);

    field = mxGetField(pars_, 0, "numlbfgs");
    if( field != NULL && mxIsDouble(field) == 1 && mxIsComplex(field) == 0 && mxGetM(field)*mxGetN(field) == 1 )
      numbfgsvecs = (size_t)mxGetScalar(field);

    field = mxGetField(pars_, 0, "auxtn");
    if( field != NULL && mxIsDouble(field) == 1 && mxIsComplex(field) == 0 && mxGetM(field)*mxGetN(field) == 1 )
      doAR = (size_t)mxGetScalar(field);

    field = mxGetField(pars_, 0, "ranktol");
    if( field != NULL && mxIsDouble(field) == 1 && mxIsComplex(field) == 0 && mxGetM(field)*mxGetN(field) == 1 )
      rankredtol = mxGetScalar(field);

    field = mxGetField(pars_, 0, "precond");
    if( field != NULL && mxIsDouble(field) == 1 && mxIsComplex(field) == 0 && mxGetM(field)*mxGetN(field) == 1 )
      precond = (size_t)mxGetScalar(field);

    field = mxGetField(pars_, 0, "fillden");
    if( field != NULL && mxIsDouble(field) == 1 && mxIsComplex(field) == 0 && mxGetM(field)*mxGetN(field) == 1 )
      gdens = mxGetScalar(field);

    field = mxGetField(pars_, 0, "reorder");
    if( field != NULL && mxIsDouble(field) == 1 && mxIsComplex(field) == 0 && mxGetM(field)*mxGetN(field) == 1 )
      reorder = (size_t)mxGetScalar(field);

    field = mxGetField(pars_, 0, "gaptol");
    if( field != NULL && mxIsDouble(field) == 1 && mxIsComplex(field) == 0 && mxGetM(field)*mxGetN(field) == 1 )
      gaptol = mxGetScalar(field);

    field = mxGetField(pars_, 0, "checkbd");
    if( field != NULL && mxIsDouble(field) == 1 && mxIsComplex(field) == 0 && mxGetM(field)*mxGetN(field) == 1 )
      checkbd = (ptrdiff_t)mxGetScalar(field);

    field = mxGetField(pars_, 0, "bdtype");
    if( field != NULL && mxIsDouble(field) == 1 && mxIsComplex(field) == 0 && mxGetM(field)*mxGetN(field) == 1 )
      typebd = (size_t)mxGetScalar(field);

    field = mxGetField(pars_, 0, "writeraw");
    if( field != NULL && mxIsDouble(field) == 1 && mxIsComplex(field) == 0 && mxGetM(field)*mxGetN(field) == 1 )
      writeraw = (size_t)mxGetScalar(field);
    else writeraw = 0;

    field = mxGetField(pars_, 0, "soln_factored");
    if( field != NULL && mxIsDouble(field) == 1 && mxIsComplex(field) == 0 && mxGetM(field)*mxGetN(field) == 1 )
      soln_factored = (size_t)mxGetScalar(field);
    else soln_factored = 0;

    field = mxGetField(pars_, 0, "forcerank");
    if( field != NULL && mxIsDouble(field) == 1 && mxIsComplex(field) == 0 &&
        mymax(mxGetM(field),mxGetN(field)) == numblk &&
        mymin(mxGetM(field),mxGetN(field)) == 1 )
      forcerank = mxGetPr(field);
    else forcerank = NULL;

    field = mxGetField(pars_, 0, "seed");
    if( field != NULL && mxIsDouble(field) == 1 && mxIsComplex(field) == 0 && mxGetM(field)*mxGetN(field) == 1 )
      seed = (size_t)mxGetScalar(field);
    /* else seed has already been set */

  }

  /* Rudimentary checking of params */

  /* Single checks */

  if(rho_f <= 0) {
    mexErrMsgTxt("sdplr: Parameter 'feastol' must be positive.");
  }

  if(rho_c <= 0) {
    mexErrMsgTxt("sdplr: Parameter 'centol' must be positive.");
  }

/*
  if(method != 1 && method != 2 && method != 3) {
    mexErrMsgTxt("sdplr: Parameter 'dir' must be 1, 2, or 3.");
  }
  */

  if(sigmafac <= 1.0) {
    mexErrMsgTxt("sdplr: Parameter 'penfac' must be greater than 1.0.");
  }

  if(rankreduce != 0 && rankreduce != 1) {
    mexErrMsgTxt("sdplr: Parameter 'reduce' must be 0 or 1.");
  }

  if(timelim <= 0) {
    mexErrMsgTxt("sdplr: Parameter 'limit' must be positive.");
  }

  if(printlevel != 0 && printlevel != 1) {
    mexErrMsgTxt("sdplr: Parameter 'printlevel' must be 0 or 1.");
  }

  if(dthresh_dim < 0) {
    mexErrMsgTxt("sdplr: Parameter 'dthresh_dim' must be nonnegative.");
  }

  if(dthresh_dens < 0.0 || dthresh_dens > 1.0) {
    mexErrMsgTxt("sdplr: Parameter 'dthres_dens' must be in [0,1].");
  }

  if(numbfgsvecs <= 0) {
    mexErrMsgTxt("sdplr: Parameter 'numlbfgs' must be a positive size_teger.");
  }

/*  if(precond != 0  &&
      precond != 1  && precond != 2  && precond !=  3 &&
      precond != 11 && precond != 12 && precond != 13 &&
      precond != 21 && precond != 22 && precond != 23 && precond != 100 && precond != 200) {
    mexErrMsgTxt("sdplr: Parameter 'precond' has an illegal value.");
  }
  */

/*
 * if(gdens < 0.0 || gdens > 1.0) {
    mexErrMsgTxt("sdplr: Parameter 'fillden' must be in [0,1].");
  }

  */

/*  if(reorder != 0 && reorder != 1) {
    mexErrMsgTxt("sdplr: Parameter 'reorder' must be 0 or 1.");
  }
  */

  if(gaptol <= 0) {
    mexErrMsgTxt("sdplr: Parameter 'gaptol' must be positive.");
  }

  if(checkbd != -1) {
    mexErrMsgTxt("sdplr: At this time, parameter 'checkbd' must be -1.");
  }

  if(checkbd != -1 && checkbd != 0 && checkbd != 1) {
    mexErrMsgTxt("sdplr: Parameter 'checkbd' must be -1, 0, or 1.");
  }

  if(typebd != 1) {
    mexErrMsgTxt("sdplr: Currently, parameter 'typebd' must equal 1.");
  }

  if(rankredtol <= 0.0) {
    mexErrMsgTxt("sdplr: Parameter 'ranktol' must be positive.");
  }

/*  if(doAR != 0 && doAR != 1) {
    mexErrMsgTxt("sdplr: Parameter 'auxtn' must be 0 or 1.");
  }

  */

  /* Pairwise checks */

/*
  if(doAR == 0 && (precond == 1 || precond == 11)) {
    mexErrMsgTxt("sdplr: Incompatible 'auxtn' and 'precond' parameters.");
  }
  */

  /* Write problem out */
  if(writeraw == 1)
    writedata_raw("temp.raw", m, numblk, blksz,
                  blktype, b, CAent,
                  CArow, CAcol, CAinfo_entptr,
                  CAinfo_rowcolptr, CAinfo_type,
                  CAinfo_storage);


  /* Get ready to run problem */

  MYCALLOC(maxranks, size_t, numblk);
  MYCALLOC(ranks, size_t, numblk);
//   printf("CAinfo_entptr[0] = %d\n", CAinfo_entptr[0]);
  getstorage(m, numblk, blksz, blktype, CAinfo_entptr, &n, &nr, maxranks);


//   printf("CAinfo_entptr[0] = %d\n", CAinfo_entptr[0]);

  /* Right now this is a hack */
  if(forcerank != NULL) {
    nr = 0;
    for(k = 0; k < numblk; k++) {
      if(blktype[k] == 's') j = blksz[k];
      else if(blktype[k] == 'd') j = 1;
      if((size_t)forcerank[k] >= 1 && (size_t)forcerank[k] <= j)
        maxranks[k] = (size_t)forcerank[k];
      nr += blksz[k]*maxranks[k];
    }
  }

  for(k = 0; k < numblk; k++) ranks[k] = maxranks[k];

  MYCALLOC(R, double, nr);
  MYCALLOC(lambda, double, m);

  srand(seed);
  for(h = 0; h < nr; h++) R[h] = (double)rand()/RAND_MAX - (double)rand()/RAND_MAX;
  pieces[0] = (double)0;
  pieces[1] = (double)0;
  pieces[2] = (double)0;
  pieces[3] = (double)0;
  pieces[4] = (double)0;
  pieces[5] = (double)0.0;
  pieces[6] = (double)1.0/n;
  pieces[7] = (double)1.0;

  /* Initialize lambda according to user's y0 */
  if(!y0_empty) {
    mat = mxGetPr(y0_);
    for(i = 0; i < m; i++) lambda[i] = mat[i];
  }

  if(!info0_empty) {

    field = mxGetField(info0_, 0, "majoriter");
    if(field != NULL) pieces[0] = (size_t)mxGetScalar(field);

    field = mxGetField(info0_, 0, "minoriter"); 
    if(field != NULL) pieces[1] = (size_t)mxGetScalar(field);

    field = mxGetField(info0_, 0, "cgs"); 
    if(field != NULL) pieces[3] = (size_t)mxGetScalar(field);

    field = mxGetField(info0_, 0, "time"); 
    if(field != NULL) pieces[5] = (double)mxGetScalar(field);

    field = mxGetField(info0_, 0, "penalty"); 
    if(field != NULL) pieces[6] = (double)mxGetScalar(field);

  }

  if(!r0_empty) {

    base = 0;
    for(k = 0; k < numblk; k++) {

      cell = mxGetCell(r0_,k);
      if(mxGetM(cell) != blksz[k]) {
        mexErrMsgTxt("sdplr: Cell of input 10 (perhaps converted from input 7) has incorrect dimension.\n");
        return;
      }
      if(mxGetN(cell) >= ranks[k]) {
        mat = mxGetPr(cell);
        for(i = base; i < base + blksz[k]*ranks[k]; i++) R[i] = mat[i-base];
      }
      else {
        mat = mxGetPr(cell);
        for(i = base; i < base + blksz[k]*mxGetN(cell); i++) R[i] = mat[i-base];
        for(i = base + blksz[k]*mxGetN(cell); i < base + blksz[k]*ranks[k]; i++) {
          R[i] = (double)rand()/RAND_MAX - (double)rand()/RAND_MAX ;
          if(R[i] > 0.0) R[i] =  R[i]*R[i];
          else           R[i] = -R[i]*R[i];
        }
      }
      base += blksz[k]*ranks[k];

    }

  }

//   printf("CAinfo_entptr[0] = %d\n", CAinfo_entptr[0]);
//    writedata_sdpa("temp.dat-s", m, numblk, blksz, blktype, b, CAent, CArow, CAcol, CAinfo_entptr, CAinfo_type);
   /* exit(0); */

  /* Run problem */
  sdplrlib(m, numblk, blksz, blktype, b, CAent, CArow, CAcol,
           CAinfo_entptr, CAinfo_type,
           numbfgsvecs, rho_f, rho_c, sigmafac, rankreduce, gaptol,
           checkbd, typebd, dthresh_dim, dthresh_dens, timelim,
           rankredtol, printlevel, R-1, lambda-1, maxranks, ranks,
           pieces); 

  /* Get ouput */

  if(nlhs > 0) {


    plhs[0] = mxCreateDoubleMatrix(n_, 1, mxREAL);
    mat = mxGetPr(plhs[0]);
    base = 0; base1 = 0;
    for(k = 0; k < numblk; k++) {

      if(soln_factored == 1) {

        /* Get R */
        plhs[0] = mxCreateCellMatrix(numblk, 1);
        base = 0;
        for(k = 0; k < numblk; k++) {
          cell = mxCreateDoubleMatrix(blksz[k],ranks[k],mxREAL);
          mat = mxGetPr(cell);
          for(i = base; i < base + blksz[k]*ranks[k]; i++) mat[i-base] = R[i];
          mxSetCell(plhs[0], k, cell);
          base += blksz[k]*ranks[k];
        }

      }
      else {

        if(blktype[k] == 'd') {
          for(i = 0; i < blksz[k]; i++)
            mat[base1 + i] = R[base + i]*R[base + i];
          base1 += blksz[k];
        }
        else if(blktype[k] == 's') {
          dsyrk_(&uplo, &trans, &(blksz[k]), &(ranks[k]), &one, R + base, &(blksz[k]), &zero, mat + base1, &(blksz[k]));
          for(i = 0; i < blksz[k]; i++)
            for(j = i+1; j < blksz[k]; j++)
              mat[base1 + j*blksz[k] + i] = mat[base1 + i*blksz[k] + j];
          base1 += blksz[k]*blksz[k];
        }
        base += blksz[k]*ranks[k];

      }

    }

    if(nlhs > 1) {

      /* Get dual variable */
      plhs[1] = mxCreateDoubleMatrix(m,1,mxREAL);
      mat = mxGetPr(plhs[1]);
      for(i = 0; i < m; i++) mat[i] = lambda[i];

      if(nlhs > 2) {

        MYCALLOC(infonames, char*, 6);
        for(i = 0; i < 6; i++)
          MYCALLOC(infonames[i], char, 9);
        strcpy(infonames[0], "majoriter");
        strcpy(infonames[1], "minoriter");
        strcpy(infonames[2], "cgs");
        strcpy(infonames[3], "time");
        strcpy(infonames[4], "penalty");

        plhs[2] = mxCreateStructMatrix(1,1,6,(const char **)infonames);
        if(plhs[2] == NULL) {
          mexErrMsgTxt("Internal error (8).");
          return;
        }

        for(i = 0; i < 6; i++)
          MYFREE(infonames[i]);
        MYFREE(infonames);

        /* Pass majoriter */
        /* field = mxCreateDoubleScalar(pieces[0]); */
        field = mxCreateDoubleMatrix(1, 1, mxREAL);
        *mxGetPr(field) = pieces[0];
        mxSetFieldByNumber(plhs[2],0,0,field);
 
        /* Pass minoriter */
        /* field = mxCreateDoubleScalar(pieces[1]); */
        field = mxCreateDoubleMatrix(1, 1, mxREAL);
        *mxGetPr(field) = pieces[1];
        mxSetFieldByNumber(plhs[2],0,1,field);

        /* Pass cgs */
        /* field = mxCreateDoubleScalar(pieces[3]); */
        field = mxCreateDoubleMatrix(1, 1, mxREAL);
        *mxGetPr(field) = pieces[3];
        mxSetFieldByNumber(plhs[2],0,2,field);

        /* Pass time */
        /* field = mxCreateDoubleScalar(pieces[5]); */
        field = mxCreateDoubleMatrix(1, 1, mxREAL);
        *mxGetPr(field) = pieces[5];
        mxSetFieldByNumber(plhs[2],0,3,field);

        /* Pass sigma */
        /* field = mxCreateDoubleScalar(pieces[6]); */
        field = mxCreateDoubleMatrix(1, 1, mxREAL);
        *mxGetPr(field) = pieces[6];
        mxSetFieldByNumber(plhs[2],0,4,field);

        /* Pass sc */
        /* field = mxCreateDoubleScalar(pieces[7]); */
        field = mxCreateDoubleMatrix(1, 1, mxREAL);
        *mxGetPr(field) = pieces[7];
        mxSetFieldByNumber(plhs[2],0,5,field);

        if(nlhs > 3) {

          if(soln_factored == 0) {

            /* Get R */
            plhs[3] = mxCreateCellMatrix(numblk, 1);
            base = 0;
            for(k = 0; k < numblk; k++) {
              cell = mxCreateDoubleMatrix(blksz[k],ranks[k],mxREAL);
              mat = mxGetPr(cell);
              for(i = base; i < base + blksz[k]*ranks[k]; i++) mat[i-base] = R[i];
              mxSetCell(plhs[3], k, cell);
              base += blksz[k]*ranks[k];
            }

          }
          else {

            plhs[3] = mxCreateString("sdplr: Factored form of solution returned in first output argument.");

          }

        }

      }
    }


  }

  /* Free memory */

  MYFREE(lambda);
  MYFREE(R);
  MYFREE(ranks);
  MYFREE(maxranks);

  MYFREE(CAcol);
  MYFREE(CArow);
  MYFREE(CAent);

  MYFREE(CAinfo_storage);
  MYFREE(CAinfo_type);
  MYFREE(CAinfo_rowcolptr);
  MYFREE(CAinfo_entptr);

  MYFREE(b);
  MYFREE(blktype);
  MYFREE(blksz);

  MYFREE(inblkptr);
  MYFREE(inblk);

  return;
}

size_t classifyApq(size_t* blksz, char* blktype, size_t* inblk, size_t* inblkptr, size_t p, size_t q, size_t A_tr, size_t* cons, size_t* blk, size_t* i, size_t* j)
{
  size_t r, c;

  if(A_tr) {
    r = q;
    c = p;
  }
  else {
    r = p;
    c = q;
  }

  *cons = r+1;
  *blk = inblk[c];
  if(blktype[*blk - 1] == 'd') {
    *i = *j = c - inblkptr[*blk] + 1;
  }
  else if(blktype[*blk - 1] == 's') {
    *i = (c - inblkptr[*blk]) % blksz[*blk - 1] + 1;
    *j = (c - inblkptr[*blk]) / blksz[*blk - 1] + 1;
  }

  return 0;
}

