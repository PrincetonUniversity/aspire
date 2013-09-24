#include "myinclude.h"
#ifdef __MEX
#include "mex.h"
#endif

size_t print_notes(void)
{
  printf("\n");
  printf("Note: Readdata_sdplr() assumes only one block!\n");
  printf("Note: Dual bounds currently only partially enabled!\n");
  printf("Note: Currently, ARPACK tol = 1/n. Appropriate?\n");
  printf("Note: ARPACK does not print error when it fails to converge.\n");
  printf("Note: No PC'ing for diagonal blocks.\n");
  printf("Note: Scaling doesn't work b/c of low-rank data matrices.\n");
  printf("\n");

  return 0;
}

size_t printparams(problemdata* data)
{
  printf("rho_f        = %.1e\n", data->rho_f);
  printf("rho_c        = %.1e\n", data->rho_c);
//   printf("method       = %d\n", data->method);
  printf("sigmafac     = %.1f\n", data->sigmafac);
  printf("rankreduce   = %d\n", data->rankreduce);
  printf("timelim      = %d\n", data->timelim);
  printf("dthresh_dim  = %d\n", data->dthresh_dim);
  printf("dthresh_dens = %.2f\n", data->dthresh_dens);

  printf("numbfgsvecs  = %d\n", data->numbfgsvecs);
//   printf("mixedtol     = %.1e\n", data->mixedtol);
//   printf("doAR         = %d\n", data->doAR);
  printf("rankredtol   = %.16e\n", data->rankredtol);

//   printf("pc           = %d\n", data->pc);
//   printf("gdens        = %.1f\n", data->gdens);
//   printf("reorder      = %d\n", data->reorder);
  printf("gaptol       = %.1e\n", data->gaptol);
  printf("checkbd      = %d\n", data->checkbd);
  printf("typebd       = %d\n", data->typebd);

  return 0;
}


size_t print_dimacs_errors(problemdata* data, double* R)
{
  int problem;
  size_t i, j, k, one=1;
  double CdotRRt, bty;
  double Cnorm, bnorm;
  double tempval;
  double *X, *S;
  size_t *colptr, *rowind;
  double d1,d2,d3,d4,d5,d6;

  d2 = d3 = 0.0;

//   printf("DIMACS error measures: ");

  // This routine does not work if C is a low-rank
  // data matrix. (Have I fixed this?)

  Aoper(data,R,R,data->X,1,1,data->vio);
  EASYDAXPY(data->m, -1.0, data->b, data->vio);
  CdotRRt = data->vio[0];

  bty = EASYDDOT(data->m, data->b, data->lambda);
  Cnorm = C_normdatamat(data);
  bnorm = fabs(data->b[ idamax_(&(data->m), data->b + 1, &one) ]);

  d1 = EASYDNRM2(data->m, data->vio)/(1.0 + bnorm);

#ifdef __ARPACK
  // Begin: Make sure S is based on y only
  for(k = 1; k <= data->m; k++)
    data->y[k] = -data->lambda[k];
  AToper(data, data->y, data->S, 1);
  // End: Make sure S is based on y only
  problem = Smineval(data,&tempval);
  d4 = mymax(0.0, -tempval/(1.0 + Cnorm));
#else
  problem = Smineval(data,&tempval);
  if(fabs(tempval + 1.0e10) > DBL_EPSILON) {
    d4 = mymax(0.0, -tempval/(1.0 + Cnorm));
  }
  else { d4 = -1.0e10; }
#endif

  d5 = (CdotRRt - bty)/(1.0 + fabs(CdotRRt) + fabs(bty));

  tempval = 2.0*EASYDDOT(data->XS_blkptr[data->numblk+1]-1, data->X, data->S);
  for(k = 1; k <= data->numblk; k++) {
    X = data->X + data->XS_blkptr[k] - 1;
    S = data->S + data->XS_blkptr[k] - 1;
    if(data->blktype[k] == SDPBLK && data->XS_blksto[k] == SPARSE) {
      colptr = data->XS_colptr[k];
      rowind = data->XS_rowind[k];
      for(j = 1; j <= data->blksz[k]; j++)
        for(i = colptr[j]; i <= colptr[j+1]-1; i++) // This is a bit exorbitant.
          if(rowind[i] == j)
            tempval -= X[i]*S[i];
    }
    else if(data->blktype[k] == SDPBLK && data->XS_blksto[k] == DENSE) {
      for(j = 1; j <= data->blksz[k]; j++)
        tempval -= X[SMATIND(j,j,data->blksz[k])]*S[SMATIND(j,j,data->blksz[k])];
    }
    else if(data->blktype[k] == DIAGBLK) {
      for(j = 1; j <= data->blksz[k]; j++)
        tempval -= X[j]*S[j];
    }
  }

  d6 = tempval/(1.0 + fabs(CdotRRt) + fabs(bty));

  if(fabs(d4 + 1.0e10) > DBL_EPSILON) printf("DIMACS error measures: %.2e %.2e %.2e %.2e %.2e %.2e\n", d1, d2, d3, d4, d5, d6);
  else printf("DIMACS error measures: %.2e %.2e %.2e (n/a) %.2e %.2e\n", d1, d2, d3, d5, d6);
  
  if(problem == -1)
    printf("Warning (ARPACK): Eigenvalue calculation failed to converge. Best estimate returned.\n");

  printf("\n"); 

  return 0;
}


size_t printheading(size_t start)
{
  if(start == 1) {
#ifdef __MEX
#ifdef __WIN32
    printf("\n");
    printf("        ***   SDPLR %s   ***\n\n", VERSION);
    printf("===========================================\n");
    printf(" major   minor         val         infeas  \n");
    printf("-------------------------------------------\n");
#else
    printf("\n");
    printf("       ***   SDPLR %s   ***\n\n", VERSION);
    printf("=========================================\n");
    printf(" major   minor        val        infeas\n");
    printf("-----------------------------------------\n");
#endif
#else
#ifdef __WIN32
    printf("\n");
    printf("             ***   SDPLR %s   ***\n\n", VERSION);
    printf("=====================================================\n");
    printf(" major   minor         val         infeas      time  \n");
    printf("-----------------------------------------------------\n");
#else
    printf("\n");
    printf("            ***   SDPLR %s   ***\n\n", VERSION);
    printf("===================================================\n");
    printf(" major   minor        val        infeas      time  \n");
    printf("---------------------------------------------------\n");
#endif
#endif
  }

  if(start == 0) {
#ifdef __MEX
#ifdef __WIN32
    printf("===========================================\n\n");
#else
    printf("=========================================\n\n");
#endif
#else
#ifdef __WIN32
    printf("=====================================================\n\n");
#else
    printf("===================================================\n\n");
#endif
#endif
  }

  fflush(stdout);

  return 1;
}



  

