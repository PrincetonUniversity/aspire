function f=slow_nufft_1d(alpha,omega,m)
%
% Compute directly the sums
%            n 
%    f(j) = sum alpha(k)*exp(2*pi*i*j*omega(k)/m)
%           k=1
% for j=-m/2,...,m/2-1 (or, if m is odd, for j=(-(m-1)/2,...,(m-1)/2)),
% where n is the length of alpha and omega.
%
% The complexity of the computation is O(nm).
%
% If m is missing then m=n;
%
% Yoel Shkolnisky, December 2006.

n=length(alpha);

if nargin<3
    m=n;
end

if (length(omega)~=n)
    error('alpha and omega should have the same length');
end

low_idx=-ceil((m-1)/2);
high_idx=floor((m-1)/2);

F = exp(2*pi*1i*[low_idx:high_idx]'*omega(:)'/m);

f = F*alpha(:);
