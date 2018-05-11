function sig_hat_kk = KN_noiseEst(ell,n,kk)
% function sig_hat_kk = KN_noiseEst(ell,n,kk)
%
% Code by Shira Kritchman and Boaz Nadler
% 2008, Weizmann Institute of Science
% --------------------------------------------
% DESCRIPTION:
% 	This function gets as input the eigenvalues of an SCM 
%   (sample covariance matrix) and an assumed value of its pseudorank,
%   and outputs an estimation of the noise variance,
%   under the assumption of uncorrelated homoscedastic noise.
%   This value is used in the algorithm for rank estimation, KN_rankEst
%
% INPUT:
%   ell     -  vector of eigenvalues of the SCM, of length p
%   n       -  number of samples
%   kk      -  assumed rank
%
% OUTPUT:
%   sig_hat_kk - Estimate of the unknown (squared) noise variance, sigma^2. 
% --------------------------------------------
% FOR MORE DETAILS SEE:
%   S. Kritchman and B. Nadler, Determining the number of components in a factor model
%   from limited noisy data, 2008
% --------------------------------------------

max_iter = 30; 
eps_threshold = 1e-5; 

p = length(ell);

sigma_0 = 1/(p-kk) * sum(ell(kk+1:end)) * 1 / (1-kk / n); 

for counter = 1:max_iter
    
    % solve quadratic equation for rho, given sigma and eigenvalues
    tmp = ell(1:kk) + sigma_0 - (p-kk)/n*sigma_0; 
    if min(tmp.^2 - 4 * ell(1:kk)*sigma_0) < 0  % otherwise get complex valued solutions
        break; 
    end
    Rho_est = zeros(kk,1); 
    Rho_est = ( tmp + sqrt( tmp.^2 - 4 * ell(1:kk)*sigma_0)  ) / 2 ; 

    if min(ell(1:kk) - Rho_est) < 0     
        fprintf('MAJOR ERROR CONSISTENT NOISE EST kk %d\n',kk); 
        ell(1:kk), Rho_est, pause; end
    Delta_l_rho = max(ell(1:kk) - Rho_est,0);
    sigma_new = 1/(p-kk) * ( sum(ell(kk+1:end)) + sum(Delta_l_rho) ) ;
    if abs(sigma_new - sigma_0)/sigma_0 < eps_threshold
        break;
    else
        sigma_0 = sigma_new; 
    end
end

% if counter > 10 
%     fprintf('iter consistent %d\n',counter); 
% end

sig_hat_kk = sigma_0;

return; 

