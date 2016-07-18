A matlab implementation of the rank estimation algorithm described in 

S. Kritchman and B. Nadler, 
Determining the number of components in a factor model from limited noisy data, Chem. Int. Lab. Sys. 2008.

To run the software use:

K = rank_est_KN (ell, n, beta) 

or 

[K sigma_hat] = rank_est_KN (ell, n, beta) 


where

ell = eigenvalues of the sample covariance matrix
n   = number of samples
beta= 1/2 for real valued / complex valued noise. 

OUTPUT: 
K   = estimate of rank
sigma_hat = estimate of noise variance. 

See help of each file for more advanced options.  