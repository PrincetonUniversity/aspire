function [mu_np sigma_np] = KN_mu_sigma(n,p,beta);
% function [mu_np sigma_np] = KN_mu_sigma(n,p,beta);
%
% Code by Shira Kritchman and Boaz Nadler
% 2008, Weizmann Institute of Science
% --------------------------------------------
% DESCRIPTION:
% 	This function computes the parameters mu_np and sigma_np
%   which are used to normalize ell_1, the largest eigenvalue of
%   a Wishart matrix. After the normalization the distribution of
%   ell_1 converges to a Tracy-Widom distribution:
%   Pr{ell_1 > (mu_np + s sigma_np)} --> F_beta(s)
%   These values are used in the algorithm for rank estimation, KN_rankEst.
%
% INPUT:
%   n       -  number of samples
%   p       -  dimension of samples
%   beta    -  indicator for real (beta=1) or complex (beta=2) valued observations
%
% OUTPUT:
%   [mu_np,sigma_np]
% --------------------------------------------
% FOR MORE DETAILS SEE:
%   S. Kritchman and B. Nadler, Determining the number of components in a factor model
%   from limited noisy data, 2008
%
% FOR MORE DETAILS ON THE COMPUTATION OF mu_np AND sigma_np SEE:
%
% I. M. Johnstone, High Dimensional Statistical Inference and Random
% Matrices, Proc. International Congress of Mathematicians, 2006.
%
% N. El Karoui, A rate of convergence result for the largest eigenvalue of
% complex white Wishart matrices, Annals of Probability, 34(6):2077-2117,
% 2006.
% --------------------------------------------

if beta==1
    mu_np = (sqrt(n-1/2) + sqrt(p-1/2))^2;
    sigma_np = sqrt(mu_np) * (1/sqrt(n-1/2) + 1 / sqrt(p-1/2) )^(1/3);
elseif beta==2
    P = min(n,p); N = max(n,p);  
    N_plus = N+1/2; P_plus = P+1/2;    
    Nm1_plus = N-1+1/2; Pm1_plus = P-1+1/2;
  
    mu_Nm1P_temp = (sqrt(Nm1_plus) + sqrt(P_plus))  ^2;
    mu_NPm1_temp = (sqrt(N_plus)   + sqrt(Pm1_plus))^2;
    sigma_Nm1P_temp = (sqrt(Nm1_plus) + sqrt(P_plus))   * (1/sqrt(Nm1_plus) + 1/sqrt(P_plus))  ^(1/3);
    sigma_NPm1_temp = (sqrt(N_plus)   + sqrt(Pm1_plus)) * (1/sqrt(N_plus)   + 1/sqrt(Pm1_plus))^(1/3);

    gamma_NP = (mu_Nm1P_temp * sigma_NPm1_temp^0.5) / ...
               (mu_NPm1_temp * sigma_Nm1P_temp^0.5);
    sigma_np = (1+gamma_NP) / (1/sigma_Nm1P_temp + gamma_NP/sigma_NPm1_temp);
    mu_np    = (1/sigma_Nm1P_temp^0.5 + ...
                1/sigma_NPm1_temp^0.5) / ...
               (1/(mu_Nm1P_temp * sigma_Nm1P_temp^0.5) + ...
                1/(mu_NPm1_temp * sigma_NPm1_temp^0.5));
end