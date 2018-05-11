function [K sigma_hat] = KN_rankEst(ell,n,beta,alpha,max_kk)
% function [K sigma_hat] = KN_rankEst(ell,n,beta,alpha,max_kk)
%
% Code by Shira Kritchman and Boaz Nadler
% 2008, Weizmann Institute of Science
% --------------------------------------------
% DESCRIPTION:
% 	This function gets as input the eigenvalues of an SCM 
%   (sample covariance matrix), and outputs an estimation of its 
%   pseudorank, under the assumption of uncorrelated homoscedastic noise.
%
% INPUT:
%   ell     -  vector of eigenvalues of the SCM, of length p
%   n       -  number of samples
%   beta    -  indicator for real (1) or complex (2) valued observations
%   alpha   -  confidence level, given in percents
%   max_kk  -  maximal possible value for pseudorank
%
% DEFAULT:
%   if nargin<4 alpha = 0.5
%   if nargin<5 max_kk = min(n,p)-1
%
% OUTPUT:
%   K - pseudorank estimation for the SCM
% sigma_hat - estimate of noise variance. 
% --------------------------------------------
% FOR MORE DETAILS SEE:
%   S. Kritchman and B. Nadler. Determining the number of components in a factor model
%   from limited noisy data, Chem. Int. Lab. Sys. 2008
% --------------------------------------------


p = length(ell);

if nargin<5, max_kk = min(n,p)-1; end
max_kk = min(max_kk,min(p,n)-1);

if nargin<4, alpha = 0.5; end
s_Wishart = KN_s_Wishart(alpha,beta);

sigma_arr = zeros(1,max_kk); 
for kk=1:max_kk  
    [mu_np sigma_np] = KN_mu_sigma(n,p-kk,beta);
    sig_hat_kk = KN_noiseEst(ell,n,kk);
    sigma_arr(kk) = sig_hat_kk; 
    at_least_kk_signals = n * ell(kk) > sig_hat_kk * (mu_np + s_Wishart * sigma_np);
    if ~at_least_kk_signals, break, end
end  % for kk=1:max_kk
K = kk-1;
if K > 0 
   sigma_hat = sigma_arr(K); 
else
   sigma_hat = sum(ell(1:p)) / p;
end
