function g=slow_nufft_t_1d(beta,x);
%
% Compute directly the sums
% The function computes the sums
%               n/2
%        g(j) = sum beta(k)*exp(i*k*x(j))
%              k=-n/2
% for j=1,...,m.
% where n is the length of beta and x.
%
% The complexity of the computation is O(nm).
%
% If m is missing then m=n;
%
% Yoel Shkolnisky, December 2006.

n=length(beta);
m=length(x);

low_idx=-ceil((n-1)/2);
high_idx=floor((n-1)/2);

G = exp(1i*x(:)*[low_idx:high_idx]);
g = G*beta;
