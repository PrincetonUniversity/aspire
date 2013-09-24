#include "myinclude.h"
#ifdef __MEX
#include "mex.h"
#endif

/* #ifdef __MEX */
/* #undef exit */
/* #define exit(er) mexErrMsgTxt("SDPLR Error\n"); */
/* #endif */

size_t copystructures(problemdata* data, size_t m, size_t numblk, size_t* blksz, char* blktype, double* b,
                   double* CAent, size_t* CArow, size_t* CAcol, size_t* CAinfo_entptr,
                   char* CAinfo_type)
{
  size_t h, i, j, k, nnz;
  char label[5] = "none\0", type;
  size_t tempsize_t, tempsize_t2;
  sparsesymmmat *Asp;
  diagmat *Adiag;
  lowrankmat *Alr;
  datamat **A;

  int lim1;
  int lim2;


  data->m = m;
  data->numblk = numblk;
  data->blksz = blksz-1;
  data->blktype = blktype-1;
  data->b = b-1;

  MYCALLOC(data->A, datamat**, data->m + 1);
  for(i = 1; i <= data->m; i++)
    MYCALLOC(data->A[i], datamat*, data->numblk + 1);
  MYCALLOC(data->C, datamat*, data->numblk + 1);

  for(j = 0; j <= data->m; j++) {

    for(i = 1; i <= data->numblk; i++) {

      nnz = CAinfo_entptr[ DATABLOCKIND(j,i,data->numblk) + 1] - CAinfo_entptr[ DATABLOCKIND(j,i,data->numblk)];
      type = CAinfo_type[ DATABLOCKIND(j,i,data->numblk) ];

      if(type == 's' || type == 'd') {
        tempsize_t = nnz;
        tempsize_t2 = -1;
      }
      else if(type == 'l') {
        tempsize_t = nnz/(data->blksz[i] + 1);
        tempsize_t2 = data->blksz[i];
      }
      else {
/* #ifdef __MEX */
/*         mexErrMsgTxt("copystructures: type not one of three possible values."); */
/* #else */
        printf("copystructures: type not one of three possible values.\n");
        exit(0);
/* #endif */
      }

      if(j == 0) A = data->C + i;
      else A = data->A[j] + i;
      createdatamat(A, type, tempsize_t, tempsize_t2, label);

    }

  }

  for(j = 0; j <= data->m; j++) {

    for(i = 1; i <= data->numblk; i++) {

      type = CAinfo_type[ DATABLOCKIND(j,i,data->numblk) ];

      if(type == 'd') {

        if(j == 0) Adiag = data->C[i]->diag;
        else Adiag = data->A[j][i]->diag;

        /* Strange use of lim1, lim2 added because for loop was being entered
         * even when lim1=0 and lim2=-1. I could only guess it was because
         * of macro expansion or something weird. */

        lim1 = CAinfo_entptr[ DATABLOCKIND(j,i,data->numblk) ]; 
        lim2 = CAinfo_entptr[ DATABLOCKIND(j,i,data->numblk) + 1] - 1;

        if(lim1 <= lim2) {

          h = 0;
          for(k = CAinfo_entptr[ DATABLOCKIND(j,i,data->numblk) ]; k <= CAinfo_entptr[ DATABLOCKIND(j,i,data->numblk) + 1] - 1; k++) {
            h++;
            if(CArow[k] != CAcol[k]) {
              printf("Error (copystructures.c): Data for diagonal matrix is not diagonal.\n");
              printf("%d %d\n", CArow[k], CAcol[k]);
              exit(0);
            }
            Adiag->ind[h] = CArow[k];
            Adiag->ent[h] = CAent[k];
          }

        }

      }
      else if(type == 's') {

        if(j == 0) Asp = data->C[i]->sp;
        else Asp = data->A[j][i]->sp;

        h = 0;
        for(k = CAinfo_entptr[ DATABLOCKIND(j,i,data->numblk) ]; (int)k <= (int)CAinfo_entptr[ DATABLOCKIND(j,i,data->numblk) + 1] - 1; k++) {
          h++;
          if(CArow[k] < CAcol[k]) { tempsize_t = CArow[k]; CArow[k] = CAcol[k]; CAcol[k] = tempsize_t; }
          Asp->row[h] = CArow[k];
          Asp->col[h] = CAcol[k];
          Asp->ent[h] = CAent[k];
        }

      }
      else if(type == 'l') {

        if(j == 0) Alr = data->C[i]->lr;
        else Alr = data->A[j][i]->lr;

        mydcopy(Alr->ncol, CAent + CAinfo_entptr[ DATABLOCKIND(j,i,data->numblk) ], 1, Alr->d + 1, 1);
        mydcopy(Alr->nrow*Alr->ncol, CAent + CAinfo_entptr[ DATABLOCKIND(j,i,data->numblk) ] + Alr->ncol, 1, Alr->ent + 1, 1);

      }

    }

  }


  /*
  for(j = 0; j <= data->m; j++) {

    if(j == 0) A = data->C;
    else A = data->A[j];

    for(i = 1; i <= data->numblk; i++) {

      if(data->blktype[i] == 's') {

        datamat_quicksort3(A[i]->sp->col, A[i]->sp->row, A[i]->sp->ent, 1, A[i]->sp->nnz);

      }
      else if(data->blktype[i] == 'd') {

        datamat_quicksort2(A[i]->diag->ind, A[i]->diag->ent, 1, A[i]->diag->nnz);

      }

    }

  }
  */


  return 0;
}

/*
size_t datamat_quicksort2(size_t* A1, double* A2, size_t p, size_t r)
{
   size_t q;

   if(p < r) 
   {
      q = datamat_partition2(A1, A2, p, r);
      datamat_quicksort2(A1, A2, p, q);
      datamat_quicksort2(A1, A2, q+1, r);
   }

   return 1;
}

size_t datamat_partition2(size_t* A1, double* A2, size_t p, size_t r)
{
   size_t i, j, t1;
   double m, k, t2;

   k = A1[p]; m = A2[p];
   i = p-1;
   j = r+1;
   while(i < j) {
      do j--;
      while(A1[j] > k || (A1[j] == k && A2[j] > m) );
      do i++;
      while(A1[i] < k || (A1[i] == k && A2[i] < m) );
      if(i < j) {
         t1 = A1[j]; t2 = A2[j];
         A1[j] = A1[i]; A2[j] = A2[i];
         A1[i] = t1; A2[i] = t2;
      }
      else return j;
   }

   return 0;
}


size_t datamat_quicksort3(size_t* A1, size_t* A2, double* A5, size_t p, size_t r)
{
   size_t q;

   if(p < r) 
   {
      q = datamat_partition3(A1, A2, A5, p, r);
      datamat_quicksort3(A1, A2, A5, p, q);
      datamat_quicksort3(A1, A2, A5, q+1, r);
   }

   return 1;
}

size_t datamat_partition3(size_t* A1, size_t* A2, double* A5, size_t p, size_t r)
{
   size_t i, j;
   size_t sv1, sv2;
   size_t t1, t2;
   double t5;

   sv1 = A1[p]; sv2 = A2[p];
   i = p-1;
   j = r+1;
   while(i < j) {
      do j--;
      while(A1[j] > sv1 || (A1[j] == sv1 && A2[j] > sv2) );
      do i++;
      while(A1[i] < sv1 || (A1[i] == sv1 && A2[i] < sv2) );
      if(i < j) {
         t1 = A1[j]; t2 = A2[j]; t5 = A5[j];
         A1[j] = A1[i]; A2[j] = A2[i]; A5[j] = A5[i];
         A1[i] = t1; A2[i] = t2; A5[i] = t5;
      }
      else return j;
   }

   return 0;
}

*/

