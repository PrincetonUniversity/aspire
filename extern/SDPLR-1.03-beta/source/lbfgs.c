#include "myinclude.h"

size_t dirlbfgs(problemdata* data, lbfgsvec* vecs, double* dir, double* grad, size_t oldest, size_t numlbfgsvecs, size_t scale)
{
  size_t i, ind;
  double beta;
  lbfgsvec *vec;

  EASYDCOPY(data->nr, grad, dir);

  for(i = 1; i <= numlbfgsvecs; i++) {
//    if(oldest - i > 0) ind = oldest - i; else ind = oldest - i + numlbfgsvecs;
    if(oldest> i) ind = oldest - i; else ind = numlbfgsvecs+oldest - i;
    vec = vecs + ind;
    vec->a = EASYDDOT(data->nr, vec->s, dir); vec->a *= vec->rho;
    EASYDAXPY(data->nr, -(vec->a), vec->y, dir);
  }

  for(i = numlbfgsvecs; i >= 1; i--) {
//    if(oldest - i > 0) ind = oldest - i; else ind = oldest - i + numlbfgsvecs;
    if(oldest> i) ind = oldest - i; else ind = numlbfgsvecs+oldest - i;
    vec = vecs + ind;
    beta = EASYDDOT(data->nr, vec->y, dir); beta *= vec->rho;
    EASYDAXPY(data->nr, vec->a - beta, vec->s, dir);
  }

  if(scale) EASYDSCAL(data->nr, -1.0, dir);
  
  return 1;
}


size_t updatelbfgs1(problemdata* data, lbfgsvec* vecs, double* grad, size_t oldest)
{
  copyscaledvectovec((vecs+oldest)->y, -1.0, grad, data->nr);
  return 1;
}

size_t updatelbfgs2(problemdata* data, lbfgsvec* vecs, double* dir, double* grad, double stepsize, size_t* oldest, size_t numlbfgsvecs)
{
  lbfgsvec *vec;

  vec = vecs + *oldest;
  
  copyscaledvectovec(vec->s, stepsize, dir, data->nr);
  EASYDAXPY(data->nr, 1.0, grad, vec->y);
    
  vec->rho = EASYDDOT(data->nr, vec->y, vec->s);
  vec->rho = 1.0/vec->rho;
  *oldest = (*oldest)%numlbfgsvecs + 1;

  return 1;
}


