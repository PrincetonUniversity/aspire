function [ K ] = MP_rankEst( D, n, var_hat )
%Use Marchenko Pastur distribution to estimate number of components
%   Input: D: eigenvalues
%          n: sample size
%          var_hat: estimated variance
%   Output: K: number of signal components.

l_D = length(D);
lambda = l_D/n;
K = length(find(D > var_hat*(1+sqrt(lambda))^2 ));

end

