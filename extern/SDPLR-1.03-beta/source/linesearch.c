#include "myinclude.h"
#ifdef __MEX
#include "mex.h"
#endif

/* #ifdef __MEX */
/* #undef exit */
/* #define exit(er) mexWarnMsgTxt("SDPLR Error\n"); return mxGetNaN(); */
/* #endif */

double linesearch(problemdata* data, double* R, double* D, double max, double* funcval, size_t update)
{
  size_t i;

  double ss=0.0;
  double quartic[5], cubic[4], roots[3], f0, fmax, f1, f2, f3, minval;

  // This code assumes that the violation for R is stored in data->vio,
  // i.e., function() has not been called since function(R) was
  // evaluated. Same for RRt.

  ZEROVEC(global_ARD, data->m); global_ARD[0] = 0.0;
  ZEROVEC(global_ADD, data->m); global_ADD[0] = 0.0;

  Aoper(data,R,D,global_UVt,0,1,global_ARD);
  EASYDSCAL(data->m,2.0,global_ARD); global_ARD[0] *= 2.0;

  Aoper(data,D,D,global_UVt,1,1,global_ADD);

  // Calculate constant for quartic
  quartic[0] = data->vio[0] -
               EASYDDOT(data->m, data->lambda, data->vio) + 
               0.5*data->sigma*pow(EASYDNRM2(data->m, data->vio), 2.0);
//   quartic[0] = *funcval;

  // Calculate linear for quartic
  quartic[1] = global_ARD[0] - EASYDDOT(data->m, data->lambda, global_ARD) +
               data->sigma*EASYDDOT(data->m, data->vio, global_ARD);

  // Calculate quadratic for quartic
  quartic[2] = global_ADD[0] - EASYDDOT(data->m, data->lambda, global_ADD) +
               data->sigma*(EASYDDOT(data->m, data->vio, global_ADD) +
               0.5*pow(EASYDNRM2(data->m, global_ARD), 2.0));

  // Calculate cubic for quartic
  quartic[3] = data->sigma*EASYDDOT(data->m, global_ARD, global_ADD);

  // Calculate quartic for quartic
  quartic[4] = 0.5*data->sigma*pow(EASYDNRM2(data->m, global_ADD), 2.0);

  // Calculate cubic
  cubic[0] = 1.0*quartic[1]; if(cubic[0] > DBL_EPSILON) {
    printf("Problem!  %f should be less then 0.\n", cubic[0]); return 0.0; }
  cubic[1] = 2.0*quartic[2];
  cubic[2] = 3.0*quartic[3];
  cubic[3] = 4.0*quartic[4];

  if(fabs(cubic[3]) < DBL_EPSILON) {
    printf("Surprise! Got a quadratic function!\n");
    exit(0);
  }

  roots[0] = roots[1] = roots[2] = 1.0e10;
  gsl_poly_solve_cubic(cubic[2]/cubic[3], cubic[1]/cubic[3], cubic[0]/cubic[3], roots, roots+1, roots+2);
  //printf("%f    %f    %f\n", roots[0], roots[1], roots[2]);

  f0 = quartic[0];

  fmax = gsl_poly_eval(quartic, 5, max);

  if(fabs(roots[0] - 1.0e10) < DBL_EPSILON ||
                   roots[0] < DBL_EPSILON ||
             roots[0] - max > DBL_EPSILON) f1 = 1.0e20;
  else f1 = gsl_poly_eval(quartic, 5, roots[0]);

  if(fabs(roots[1] - 1.0e10) < DBL_EPSILON || roots[1] < DBL_EPSILON || roots[1] - max > DBL_EPSILON) f2 = 1.0e20;
  else f2 = gsl_poly_eval(quartic, 5, roots[1]);

  if(fabs(roots[2] - 1.0e10) < DBL_EPSILON || roots[2] < DBL_EPSILON || roots[2] - max > DBL_EPSILON) f3 = 1.0e20;
  else f3 = gsl_poly_eval(quartic, 5, roots[2]);

  minval = 1.0e20;
  minval = mymin(f0,minval);
  minval = mymin(fmax,minval);
  minval = mymin(f1,minval);
  minval = mymin(f2,minval);
  minval = mymin(f3,minval);

  if(fabs(f0 - minval) < DBL_EPSILON) ss = 0.0;
  if(fabs(fmax - minval) < DBL_EPSILON) ss = max;
  if(fabs(f1 - minval) < DBL_EPSILON) ss = roots[0];
  if(fabs(f2 - minval) < DBL_EPSILON) ss = roots[1];
  if(fabs(f3 - minval) < DBL_EPSILON) ss = roots[2];

  // Setup next values

  *funcval = minval;
  if(update) {
    for(i = 1; i <= data->m; i++)
      data->vio[i] += ss*(global_ARD[i] + ss*global_ADD[i]);
    data->vio[0] += ss*(global_ARD[0] + ss*global_ADD[0]);
  }

  return ss;
}

