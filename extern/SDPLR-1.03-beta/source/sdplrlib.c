//! \page sdplrlib The Main Algorithm (sdplrlib.c)
//! 
//! This page describes the main algorithm of SDPLR.

#include "myinclude.h"
#include "globalsmain.h"
#ifdef __MEX
#include "mex.h"
#endif

/* #ifdef __MEX */
/* #undef exit */
/* #define exit(er) mexWarnMsgTxt("SDPLR Error\n"); goto END_MAJOR_ITERATIONS; */
/* #endif */

#define PRINTFREQ      60.0

#define LAMBDAUPDATECT 1
#define SIGMASTRATEGY  1 // If 1, then LAMBDAUPDATECT should equal 1.

#define BEGIN 1
#define END   0

#define EASY   'e'
#define MEDIUM 'm'
#define HARD   'h'

void myprint(size_t majiter, size_t iter, double val, double rho_f_val, size_t CG, double totaltime);
size_t do_scaling(problemdata *data, double value, double *norm);

size_t sdplrlib (size_t m, size_t numblk, ptrdiff_t *blksz, char *blktype, double *b,
              double *CAent, size_t *CArow, size_t *CAcol, size_t *CAinfo_entptr,
              char *CAinfo_type, size_t numbfgsvecs, double rho_f, double
              rho_c, double sigmafac, size_t rankreduce, double gaptol, ptrdiff_t
              checkbd, size_t typebd, size_t dthresh_dim, double dthresh_dens,
              size_t timelim, double rankredtol, size_t printlevel, double *R,
              double *lambda, size_t* maxranks, size_t *ranks, double *pieces)
{
  // Paramters that are stored in 'pieces'
  size_t    majiter, iter, CG, curr_CG, lambdaupdate;
  double timeoffset;

  // Algorithm declarations
  size_t          lambdaupdatect, oldest;
  double       val, rho_c_val, rho_f_val;
  double       alpha, rho_c_tol, normb, normC;
  double      *D;
  double       bestbd;
  lbfgsvec    *vecs=NULL;
  problemdata *data;

  // Various misc declarations
  size_t  i, k, localiter=0, recalc=5, recalcct=10000000;
  char line[1000], difficulty='h';
  double bestinfeas=1.0e10;
  double sc, overallsc, tv; // scaling related
  double lastval; // monitor slowing down of optimization
  double origval;

  // Timing declarations and initialization
  double timeprintfreq = 0.0;
#ifdef __WIN32
  clock_t             timeorig;
  timeorig = clock();
#else
  double              timeorig;
  struct tms          timearray;
  times(&timearray);
  timeorig = timearray.tms_utime;
#endif
  timeoffset = 0.0;

  // Allocate space for main data
  MYCALLOC(data, problemdata, 1);

  // Copy input size_to data
  copystructures(data, m, numblk, blksz, blktype, b, CAent, CArow, CAcol,
                 CAinfo_entptr, CAinfo_type);

  // Set passed variables to data and global
  data->lambda = lambda;

  // Copy user params size_to data
  data->numbfgsvecs  = numbfgsvecs;
  data->rho_f        = rho_f;
  data->rho_c        = rho_c;
  data->sigmafac     = sigmafac;
  data->rankreduce   = rankreduce;
  data->gaptol       = gaptol;
  data->checkbd      = checkbd;
  data->typebd       = typebd;
  data->dthresh_dim  = dthresh_dim;
  data->dthresh_dens = dthresh_dens;
  data->timelim      = timelim;
  data->rankredtol   = rankredtol;
  data->printlevel   = printlevel;

  // Initialize all the data structures
  initialize(data, maxranks);


  /*** BEGIN: Create data structures needed in this function. (Note: Others created in initialize) ***/

  // Direction
  MYCALLOC (D, double, data->nr + 1);

  // Exact linesearch
  MYCALLOC(global_UVt, double, data->XS_blkptr[data->numblk+1]-1 + 1);
  MYCALLOC(global_ARD, double, data->m + 1);
  MYCALLOC(global_ADD, double, data->m + 1);

  // LBFGS
  MYCALLOC (vecs, lbfgsvec, data->numbfgsvecs + 1);
  for (i = 1; i <= data->numbfgsvecs; i++) {
    MYCALLOC ((vecs + i)->s, double, data->nr + 1);
    MYCALLOC ((vecs + i)->y, double, data->nr + 1);
  }

  /*** END: Create data structures needed in this function. (Note: Others created in initialize) ***/


  // Setup algorithm parameters, quantities
  i = 1; normb          = fabs (data->b[idamax_ (&(data->m), data->b + 1, &i)]);
  normC                 = C_normdatamat (data);
  lambdaupdatect        = LAMBDAUPDATECT;
  oldest                = 1;
  bestbd                = -1.0e+20;

  /*** Now really get started ***/
  /*** Now really get started ***/

  // printparams(data);
  if(data->printlevel > 0) printheading(BEGIN);

  // Setup ranks as read in from calling routine
  data->nr = 0;
  for (k = 1; k <= data->numblk; k++) {
    data->rank[k] = ranks[k - 1];
    data->nr += data->blksz[k] * data->rank[k];
  }

  // Get algorithm status parameters from calling routine
  majiter      = (size_t)    pieces[0];
  iter         = (size_t)    pieces[1];
  lambdaupdate = (size_t)    pieces[2];
  CG           = (size_t)    pieces[3];
  curr_CG      = (size_t)    pieces[4];
  timeoffset   = (double) pieces[5];
  data->sigma  = (double) pieces[6];
  overallsc    = (double) pieces[7];

  // Do scaling and inititalize total time correctly
  if(SCALE_OBJ && normC - 1.0e-10 > DBL_EPSILON) do_scaling(data,overallsc,&normC);
//   printf("overallsc = %f\n", overallsc);
//   printf("normC = %f\n", normC);
  data->totaltime = timeoffset;
  
  // Now can setup rho_c_tol
  rho_c_tol = data->rho_c / data->sigma;

  // Calculate val, grad, etc.
  essential_calcs (data, R, normC, normb, &val, &rho_c_val, &rho_f_val);
//   printf("val = %f\n", val);
  
  // Save original value
  origval = overallsc*val;

#ifdef __MEX
        printf(""); fflush(stdout);
        mexEvalString("drawnow;");
#endif

//  printf("Hey 0\n"); fflush(stdout); mexEvalString("drawnow;");

  /*** ITERATE ITERATE ITERATE ***/
  /*** ITERATE ITERATE ITERATE ***/
  /*** ITERATE ITERATE ITERATE ***/

  while (majiter++ < 100000) {

    while ( (!SIGMASTRATEGY && lambdaupdate < lambdaupdatect) || (SIGMASTRATEGY && difficulty != EASY) ) {

      // Increment lambda counter, reset local iter count and lastval
      lambdaupdate++; localiter = 0; lastval = 1.0e10;

      // Break if already meeting stopping criterion
      if (rho_c_val <= rho_c_tol) break;

      while(rho_c_val - rho_c_tol > DBL_EPSILON) { 

        // Increment both iter and localiter counts
        iter++; localiter++;

        // Direction calculation
        copyscaledvectovec (D, -1.0, data->G, data->nr);
        dirlbfgs(data, vecs, D, data->G, oldest, data->numbfgsvecs, 1);
        updatelbfgs1(data, vecs, data->G, oldest);

        // If somehow we don't have descent, revert to steepest descent
        if (EASYDDOT (data->nr, D, data->G) >= 0.0)
          copyscaledvectovec (D, -1.0, data->G, data->nr);

        // Linesearch plus variable update
        lastval = val;
        alpha = linesearch (data, R, D, 1.0, &val, 1);
        EASYDAXPY (data->nr, alpha, D, R);

        // Refresh all the essentials
        if(recalc == 0) {
          essential_calcs (data, R, normC, normb, &val, &rho_c_val, &rho_f_val);
          recalc = recalcct;
        }
        else {
          gradient(data, R);
          rho_c_val = EASYDNRM2(data->nr, data->G)/(1.0 + normC);
          rho_f_val = EASYDNRM2(data->m, data->vio)/(1.0 + normb);
          recalc--;
        }

        // Direction calculation (necessary post-processing)
        if(data->numbfgsvecs > 0)
          updatelbfgs2 (data, vecs, D, data->G, alpha, &oldest, data->numbfgsvecs);

        // If PRINTFREQ seconds have passed since last major iteration, print an update
        timeprintfreq += current_time(timeorig) + timeoffset - data->totaltime;
        if(timeprintfreq - PRINTFREQ > DBL_EPSILON) {
          myprint(-1, iter, overallsc*val, rho_f_val, CG, data->totaltime);
          timeprintfreq -= PRINTFREQ;
        }

        // Update totaltime
//         printf("current_time = %f\n", current_time(timeorig));
//         printf("timeoffset = %f\n", timeoffset);
        data->totaltime = current_time(timeorig) + timeoffset;
//         printf("current_time = %f\n", current_time(timeorig));
//         printf("current_time = %f\n", current_time(timeorig));
//         printf("current_time = %f\n", current_time(timeorig));
//         printf("timeoffset = %f\n", timeoffset);
//         printf("data->totaltime = %f   %f\n", data->totaltime, (double)current_time(timeorig) + (double)timeoffset);

        // Possibly terminate
        if (data->totaltime >= data->timelim || rho_f_val <= data->rho_f || iter >= 10000000 || CG >= 10000000) {
            EASYDAXPY (data->m, -data->sigma, data->vio, data->lambda);
            essential_calcs (data, R, normC, normb, &val, &rho_c_val, &rho_f_val);
            goto END_CURR_MAJOR_ITERATION;
          }
        
        //printf("iter = %d  CGS = %4d  val = %f   ss = %.1e   \n", iter, CG, val, alpha);

        bestinfeas = mymin(rho_f_val, bestinfeas);

      }

#ifdef __MEX
        printf(""); fflush(stdout);
        mexEvalString("drawnow;");
#endif


      // Update Lagrange multipliers and recalculate essentials
      EASYDAXPY (data->m, -data->sigma, data->vio, data->lambda);

      tv = EASYDNRM2(data->m,data->lambda);
      if(SCALE_OBJ && normC - 1.0e-10 > DBL_EPSILON && majiter >= 2 && (tv - 10.0 > DBL_EPSILON || DBL_EPSILON < 0.1 - tv)) {
        if(tv - 10.0 > DBL_EPSILON) sc = ( 1.0 - 0.9*pow(10.0/tv,0.1) )*tv;
        else                        sc = ( 9.0*pow(10.0*tv,0.1) + 1.0 )*tv;
        overallsc *= sc;
        do_scaling(data,sc,&normC);
      }

      essential_calcs (data, R, normC, normb, &val, &rho_c_val, &rho_f_val);

      // Decide how difficult this particular subproblem was
      if(SIGMASTRATEGY) {
          if(localiter <= 10)                           difficulty = EASY;
          else if(10+1 <= localiter && localiter <= 50) difficulty = MEDIUM;
          else if(50+1 <= localiter)                    difficulty = HARD;
      }

    } // end single major iteration

    // If user wants, calculate dual bound and save best so far.
    // Assumes only one block in SDP and dual shift is by identity
    // if (data->checkbd == 1) bestbd = mymax (bestbd, dualbound (data));
    
    // Check found grossly unbound value (indicative of infeasibility)
    if(overallsc*val - 1.0e10*fabs(origval) > DBL_EPSILON) {
      printf("Cannot reduce infeasibility any further.\n");
      goto END_MAJOR_ITERATIONS;
    }

  END_CURR_MAJOR_ITERATION:

    if(data->printlevel > 0) {
      if (data->checkbd == 1) {
        sprintf (line, "%3d %6d % .7e %.1e % .7e %5d    [ %5d  %5d ]\n", majiter, iter, val, rho_f_val, bestbd, (size_t) data->totaltime, curr_CG, CG);
        printf ("%s",line); fflush (stdout);
      }
      if (data->checkbd == 0 || data->checkbd == -1) myprint(majiter, iter, overallsc*val, rho_f_val, CG, data->totaltime);
    }

#ifdef __MEX
    mexEvalString("drawnow;");
#endif

    if (data->totaltime >= data->timelim || rho_f_val <= data->rho_f || iter >= 10000000 || CG >= 10000000)
      goto END_MAJOR_ITERATIONS;

    if (_isnan (val)) { printf ("Error(sdplrlib): Got NaN.\n"); return 0; } // Sam

    /*** Get ready for next major iteration ***/

    // Adjust rank (hopefully reduce)
    if (data->rankreduce == 1) dorankreduce (data, R);

    // Update sigma
    do {
      data->sigma *= data->sigmafac;
      essential_calcs (data, R, normC, normb, &val, &rho_c_val, &rho_f_val);
      rho_c_tol = data->rho_c / data->sigma;
    } while (rho_c_val <= rho_c_tol);

    // Refresh some parameters
    lambdaupdate = 0;
    curr_CG = 0;
    if(SIGMASTRATEGY) difficulty = 'h';
    timeprintfreq = 0.0;

    // Clear BFGS vectors
    for (i = 1; i <= data->numbfgsvecs; i++) {
      ZEROVEC ((vecs + i)->s, data->nr);
      ZEROVEC ((vecs + i)->y, data->nr);
    }

  } // end major iterations


END_MAJOR_ITERATIONS:

  // Final Lanczos
//   if (data->checkbd == 0) {
//     bestbd = mymax (bestbd, dualbound (data));
//     data->totaltime = current_time (timeorig) + timeoffset;
//     printf ("%.7e  %d\n", bestbd, (size_t) data->totaltime);
//   }

  if(data->printlevel > 0) {
    printheading(END);
    essential_calcs (data, R, normC, normb, &val, &rho_c_val, &rho_f_val);
    print_dimacs_errors (data, R);
  }

  // Return info to calling program
  pieces[0] = (double) majiter;
  pieces[1] = (double) iter;
  pieces[2] = (double) lambdaupdate;
  pieces[3] = (double) CG;
  pieces[4] = (double) curr_CG;
  pieces[5] = (double) data->totaltime;
  pieces[6] = (double) data->sigma;
  pieces[7] = (double) overallsc;
  for (k = 1; k <= data->numblk; k++)
    ranks[k - 1] = data->rank[k];


  // Clean up items allocated in this file
  MYFREE (global_UVt);
  MYFREE (global_ARD);
  MYFREE (global_ADD);
  for (i = 1; i <= data->numbfgsvecs; i++) {
    MYFREE ((vecs + i)->s);
    MYFREE ((vecs + i)->y);
  }
  MYFREE (vecs);
  MYFREE (D);
  deinitialize (data);
  MYFREE (data); 

  return 0;
}

void myprint(size_t majiter, size_t iter, double val, double rho_f_val, size_t CG, double totaltime)
{
  char line[1000];
#ifndef __MEX
  if(majiter < 0) sprintf(line, "       %7d  % .8e  %.1e  %6d",          iter, val, rho_f_val, (size_t)totaltime);
  else            sprintf(line, "  %3d  %7d  % .8e  %.1e  %6d", majiter, iter, val, rho_f_val, (size_t)totaltime);
#else
  if(majiter < 0) sprintf(line, "       %7d  % .8e  %.1e",          iter, val, rho_f_val);
  else            sprintf(line, "  %3d  %7d  % .8e  %.1e", majiter, iter, val, rho_f_val);
#endif
  printf("%s",line);
  printf("\n");
  fflush(stdout);
}

size_t do_scaling(problemdata *data, double value, double *norm)
{
  size_t j, k;

  for(k = 1; k <= data->numblk; k++) {
    if(data->blktype[k] == 's') {
      if(data->C[k]->type == 's') {
        for(j = 1; j <= data->C[k]->sp->nnz; j++) {
          data->C[k]->sp->ent[j] /= value;
        }
      }
      else if(data->C[k]->type == 'l') 
        EASYDSCAL(data->C[k]->lr->ncol, 1.0/value, data->C[k]->lr->d);
    }
    else if(data->blktype[k] == 'd')
      for(j = 1; j <= data->C[k]->diag->nnz; j++) 
        data->C[k]->diag->ent[j] /= value;
  }

  for(j = data->AA_rowptr[0]; j <= data->AA_rowptr[1]-1; j++) {
    data->AA_colval_one[j] /= value;
    data->AA_colval_two[j] /= value;
  }

  *norm = C_normdatamat (data);

  EASYDSCAL(data->m, 1.0/value, data->lambda);

  return 0;

}
