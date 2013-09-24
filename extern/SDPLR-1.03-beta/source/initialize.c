#include "myinclude.h"
#ifdef __MEX
#include "mex.h"
#endif

// size_t locatetype(problemdata *data, size_t blk, char type, size_t **passed_ind, size_t *passed_nind);

size_t initialize (problemdata * data, size_t * maxranks)
{
  size_t                 h, i, j, k, ct;
  /* size_t  p, q; */
  size_t               **AGGMAT;
  size_t                *colptr;
  sparsesymmmat      *A;
  /* size_t                *xadj, *adjncy, numedges, dim, options[10]; */
  /* size_t                *perm, *invp; */
  size_t                *spind, nspind;
  /* size_t  *temprow, *tempcol; */
  /* char                type; */
  /* lowrankmat         *LR; */
  // Temp variables
  datamat* dA;
  size_t col, row, ctr;


  // The following three sections of code are very similar
  // and could be simplified.

  // Setup pointers to find lowrank data matrices. Used in several routines.
  // Right now, we do not set this up for diagonal blocks.
  MYCALLOC (data->lrind, size_t *, data->numblk + 1);
  MYCALLOC (data->nlrind, size_t, data->numblk + 1);
  for (k = 1; k <= data->numblk; k++)
    if (data->blktype[k] == 's') 
      locatetype(data, k, 'l', &(data->lrind[k]), &(data->nlrind[k]));

  // BEGIN: Setup new data structures
  // BEGIN: Setup new data structures
  // BEGIN: Setup new data structures
  // BEGIN: Setup new data structures
  
  // Setup structure of X- and S-type matrices.
  
  // For non-diag blocks, this uses somewhat exorbitant allocation of
  // an entire blksz[k] x blksz[k] lower triangular size_teger matrix,
  // but this is the easiest thing right now.
 
  // Allocate space for block related posize_ters
  
  MYCALLOC(data->XS_blkptr, size_t , data->numblk+2);
  MYCALLOC(data->XS_blksto, char, data->numblk+1);
  MYCALLOC(data->XS_colptr, size_t*, data->numblk+1);
  MYCALLOC(data->XS_rowind, size_t*, data->numblk+1);

  // First posize_ter is to position 1

  data->XS_blkptr[1] = 1;

  // Loop over blocks

  for (k = 1; k <= data->numblk; k++) {

    // First find location of sparse data matrices (i.e., not low-rank
    // or user) for this block.
   
    locatetype(data, k, 's', &spind, &nspind);

    // Condition on block type

    if(data->blktype[k] == DIAGBLK) {

      // If block type is diagonal, there is no need for colptr or
      // rowind. Just set blkptr and blksto.

      data->XS_blkptr[k+1] = data->XS_blkptr[k] + data->blksz[k];
      data->XS_blksto[k] = DENSE;

      // Determine absolute position of nonzero ntries of k-th block of
      // Ai in k-th block of S.

      for(i = 0; i <= data->m; i++) {
        if(i == 0) dA = data->C[k];
        else       dA = data->A[i][k];
        if(dA->type == 'd')
          for(j = 1; j <= dA->diag->nnz; j++)
            dA->diag->XS_in[j] = data->XS_blkptr[k] - 1 + dA->diag->ind[j];
            // printf("A(i=%d,k=%d,r=%d,c=%d) --> %d\n", i, k, dA->diag->ind[j], dA->diag->ind[j], data->XS_blkptr[k] - 1 + dA->diag->ind[j]);
      }


    }
    else if(data->blktype[k] == SDPBLK) {

      // If block type is SDP...

      // First thing is to determine whether storage will be sparse
      // or dense. So we figure out positions of aggregate nonzeros
      // (lower triangular part only) and count them.

      MYCALLOC(AGGMAT, size_t*, data->blksz[k] + 1);
      for(i = 1; i <= data->blksz[k]; i++) MYCALLOC(AGGMAT[i], size_t, i + 1);
      for(i = 1; i <= data->blksz[k]; i++) for(j = 1; j <= i; j++) AGGMAT[i][j] = 0;

      for (i = 1; i <= nspind; i++) {
        if(spind[i] == 0) A = data->C          [k]->sp;
        else              A = data->A[spind[i]][k]->sp;
        for(j = 1; j <= A->nnz; j++) AGGMAT[A->row[j]][A->col[j]] = 1;
      }

      ct = 0;
      for(j = 1; j <= data->blksz[k]; j++) for(i = j; i <= data->blksz[k]; i++) if(AGGMAT[i][j] == 1) ct++;

      // Make decision about storage type

      if ((double) 2.0 * ct / (data->blksz[k] * (data->blksz[k] + 1)) >= data->dthresh_dens || data->blksz[k] <= data->dthresh_dim) {

        // Storage is dense. Just set blkptr and blksto.
        
        data->XS_blkptr[k+1] = data->XS_blkptr[k] + data->blksz[k]*data->blksz[k];
        data->XS_blksto[k] = DENSE;

        // Determine absolute position of nonzero ntries of k-th block
        // of Ai in k-th block of S.

        for(i = 0; i <= data->m; i++) {
          if(i == 0) dA = data->C[k];
          else       dA = data->A[i][k];
          if(dA->type == 's')
            for(j = 1; j <= dA->sp->nnz; j++) {
              row = dA->sp->row[j];
              col = dA->sp->col[j];
              dA->sp->XS_in[j] = data->XS_blkptr[k] - 1 + SMATIND(row,col,data->blksz[k]); // code assumes row >= col
            }
        }
 
      }
      else {

        // Storage is sparse, so setup colptr

        MYCALLOC(data->XS_colptr[k], size_t, data->blksz[k]+2);

        // For convenience in coding, save colptr locally
        colptr = data->XS_colptr[k];

        // Next determine colptr values

        colptr[1] = 1;
        for(j = 1; j <= data->blksz[k]; j++) {
          ct = 0;
          for(i = j; i <= data->blksz[k]; i++) if(AGGMAT[i][j] == 1) ct++;
          colptr[j+1] = colptr[j] + ct;
        }

        // Now we can allocate rowind
        
        MYCALLOC(data->XS_rowind[k], size_t, colptr[data->blksz[k]+1]-1 + 1); 

        // Setup rowind
                
        for(j = 1; j <= data->blksz[k]; j++) {
          ct = 0;
          for(i = j; i <= data->blksz[k]; i++) if(AGGMAT[i][j] == 1) {
            data->XS_rowind[k][colptr[j] + ct] = i;
            ct++;
          }
        }

        // Finally, don't forget to set blkptr and blksto. Based on
        // total number of nonzeros in this block.

        data->XS_blkptr[k+1] = data->XS_blkptr[k] + colptr[data->blksz[k]+1]-1;
        data->XS_blksto[k] = SPARSE;

        // Experimental
        // Determine absolute position of nonzero ntries of k-th block of
        // Ai in k-th block of S
        
        ct = 0;
        for(j = 1; j <= data->blksz[k]; j++) for(i = j; i <= data->blksz[k]; i++) if(AGGMAT[i][j] == 1) AGGMAT[i][j] = ++ct;

        for(i = 0; i <= data->m; i++) {
          if(i == 0) dA = data->C[k];
          else       dA = data->A[i][k];
          if(dA->type == 's')
            for(j = 1; j <= dA->sp->nnz; j++) {
              row = dA->sp->row[j];
              col = dA->sp->col[j];
              dA->sp->XS_in[j] = data->XS_blkptr[k] - 1 + AGGMAT[row][col];
              // printf("A(i=%d,k=%d,r=%d,c=%d) --> %d\n", i, k, row, col, data->XS_blkptr[k] - 1 + AGGMAT[row][col]);
            }
        }

      }

      // Free AGGMAT for next round.

      for(i = 1; i <= data->blksz[k]; i++) { MYFREE(AGGMAT[i]); AGGMAT[i] = NULL; }
      MYFREE(AGGMAT); AGGMAT = NULL;

    }

    // Free location of sparse data matrices
    MYFREE(spind);

  }

  // Can now setup X and S

  MYCALLOC(data->S, double, data->XS_blkptr[data->numblk+1]-1 + 1);
  MYCALLOC(data->X, double, data->XS_blkptr[data->numblk+1]-1 + 1);

  // Also setup y

  MYCALLOC(data->y, double, data->m + 1);

  // Construct big, sparse A matrix (including C in 0-th row)
  // Watch out! A potentially weird mixture of Fortran and C indexing.
  
  // First allocate AA_rowptr

  MYCALLOC(data->AA_rowptr, size_t, data->m+2); 

  // Now determine AA_rowptr

  data->AA_rowptr[0] = 1;

  for(i = 0; i <= data->m; i++) {

    ctr = 0;

    for(k = 1; k <= data->numblk; k++) {

      if(i == 0) dA = data->C[k];
      else       dA = data->A[i][k];

           if(dA->type == 's') ctr += dA->sp  ->nnz;
      else if(dA->type == 'd') ctr += dA->diag->nnz; 

    }

    data->AA_rowptr[i+1] = data->AA_rowptr[i] + ctr;

  }

  // Next job is to allocate AA_colind and AA_colval_one and _two
  
  MYCALLOC(data->AA_colind, size_t, data->AA_rowptr[data->m+1]-1 + 1);
  MYCALLOC(data->AA_colval_one, double, data->AA_rowptr[data->m+1]-1 + 1);
  MYCALLOC(data->AA_colval_two, double, data->AA_rowptr[data->m+1]-1 + 1);

  // Then we fill them (will eventually want to sort this)

  for(i = 0; i <= data->m; i++) {

    ctr = 0;

    for(k = 1; k <= data->numblk; k++) {

      if(i == 0) dA = data->C[k];
      else       dA = data->A[i][k];

      if(dA->type == 's') {
        for(j = 1; j <= dA->sp->nnz; j++) {
          data->AA_colind    [data->AA_rowptr[i] + ctr] = dA->sp->XS_in[j];
          data->AA_colval_one[data->AA_rowptr[i] + ctr] = dA->sp->ent[j];
          if(dA->sp->row[j] == dA->sp->col[j])
            data->AA_colval_two[data->AA_rowptr[i] + ctr] =     dA->sp->ent[j];
          else
            data->AA_colval_two[data->AA_rowptr[i] + ctr] = 2.0*dA->sp->ent[j];
          ctr++;
        }
      }
      else if(dA->type == 'd') {
        for(j = 1; j <= dA->diag->nnz; j++) {
          data->AA_colind    [data->AA_rowptr[i] + ctr] = dA->diag->XS_in[j];
          data->AA_colval_one[data->AA_rowptr[i] + ctr] = dA->diag->ent[j];
          data->AA_colval_two[data->AA_rowptr[i] + ctr] = dA->diag->ent[j];
          ctr++;
        }
      }

    }

  }

  // Setup posize_ters to find lowrank data matrices. Used in several
  // routines. Right now, we do not set this up for diagonal blocks
 
  ct = 0;
  for (k = 1; k <= data->numblk; k++)
    if (data->blktype[k] == 's') {
                                     if (data->C   [k]->type == 'l') ct++;
      for (i = 1; i <= data->m; i++) if (data->A[i][k]->type == 'l') ct++;
    }
 
  data->lr_num = ct;

  MYCALLOC(data->lr_mat, size_t, data->lr_num + 1);
  MYCALLOC(data->lr_blk, size_t, data->lr_num + 1);

  ct = 0;
  for (k = 1; k <= data->numblk; k++) 
    if (data->blktype[k] == 's') {
                                     if (data->C   [k]->type == 'l') { data->lr_mat[++ct] = 0; data->lr_blk[ct] = k; }
      for (i = 1; i <= data->m; i++) if (data->A[i][k]->type == 'l') { data->lr_mat[++ct] = i; data->lr_blk[ct] = k; }
    }

  if(ct != data->lr_num) { printf("Problem getting lowrank matrices. (%d != %d)\n", ct, data->lr_num); exit(0); }

  // END: Setup new data structures
  // END: Setup new data structures
  // END: Setup new data structures
  // END: Setup new data structures

  // setup other data structures
  MYCALLOC (data->vio, double, data->m + 1);
  MYCALLOC (data->rank, size_t, data->numblk + 1);
  MYCALLOC (data->maxrank, size_t, data->numblk + 1);

  // Calculate ranks for nondiag blocks

  for (k = 1; k <= data->numblk; k++)
    data->rank[k] = data->maxrank[k] = maxranks[k-1];

  // Setup global dimension n
  data->n = 0;
  for (k = 1; k <= data->numblk; k++)
    data->n += data->blksz[k];

  // Setup global dimension nr
  data->nr = 0;
  for (k = 1; k <= data->numblk; k++)
    data->nr += data->blksz[k] * data->rank[k];

  // Now can setup data->G for gradient
  MYCALLOC (data->G, double, data->nr + 1);

  // Create global structures
  i = -1;
  for(k = 1; k <= data->numblk; k++)
    if(data->blktype[k] == SDPBLK && data->XS_blksto[k] == DENSE)
      i = mymax(i, data->blksz[k]);

  // Next setup global temp space for evaluation of function with lowrank matrices.

  ct = 0;
  for(h = 1; h <= data->lr_num; h++) {
    k = data->lr_blk[h];
    i = data->lr_mat[h];
    if(i == 0) ct = mymaxint(ct, data->rank[k]*data->C   [k]->lr->ncol);
    else       ct = mymaxint(ct, data->rank[k]*data->A[i][k]->lr->ncol);
  }

  MYCALLOC(global_UtB, double, ct+1);
  MYCALLOC(global_VtB, double, ct+1);

  return 1;
}


size_t deinitialize (problemdata * data)
{
  size_t                 i, j, k;

  MYFREE(global_UtB);
  MYFREE(global_VtB);

  MYFREE (data->vio);
  MYFREE (data->G);
  for (j = 1; j <= data->numblk; j++)
    {
      MYFREE (data->lrind[j]);
      destroydatamat (data->C[j]);
    }
  for (i = 1; i <= data->m; i++) {
    for (j = 1; j <= data->numblk; j++)
      destroydatamat (data->A[i][j]);
    MYFREE (data->A[i]);
  }

  for (k = 1; k <= data->numblk; k++) 
    if(data->blktype[k] == SDPBLK && data->XS_blksto[k] == SPARSE) {
      MYFREE(data->XS_colptr[k]);
      MYFREE(data->XS_rowind[k]); 
    }

  MYFREE(data->XS_blkptr);
  MYFREE(data->XS_blksto);
  MYFREE(data->XS_colptr);
  MYFREE(data->XS_rowind);

  MYFREE(data->S);
  MYFREE(data->X);

  MYFREE(data->y);

  MYFREE(data->AA_rowptr);
  MYFREE(data->AA_colind);
  MYFREE(data->AA_colval_one);
  MYFREE(data->AA_colval_two);

  MYFREE(data->lr_mat);
  MYFREE(data->lr_blk);

  MYFREE (data->lrind);
  MYFREE (data->nlrind);
  MYFREE (data->usind);
  MYFREE (data->nusind);
  MYFREE (data->rank);
  MYFREE (data->maxrank);
  MYFREE (data->C);
  MYFREE (data->A);

  return 1;
}


size_t locatetype(problemdata *data, size_t blk, char type, size_t **passed_ind, size_t *passed_nind)
{
  size_t i, k;
  size_t ct;
  size_t *ind, nind;

  k = blk;

  ct = 0;                        if (data->C   [k]->type == type) ct++;
  for (i = 1; i <= data->m; i++) if (data->A[i][k]->type == type) ct++;

  nind = ct;
  MYCALLOC (ind, size_t, nind + 1);

  ct = 0;                        if (data->C   [k]->type == type) ind[++ct] = 0;
  for (i = 1; i <= data->m; i++) if (data->A[i][k]->type == type) ind[++ct] = i;
  if (ct != nind) { printf ("locatetype: problem with setting up ind\n"); exit (0); }

  *passed_ind = ind;
  *passed_nind = nind;

  return 0;
}
