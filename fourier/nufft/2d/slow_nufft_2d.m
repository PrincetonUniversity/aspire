function f=slow_nufft_2d(alpha,omega,M)
%
% Compute directly the sums
%                n
%        f(j) = sum alpha(k)*exp(2*pi*i*(j1*omega1(k)+j2*omega2(k))/M)
%               k=1
% for j1,j2=-M/2,...,M/2, where n is the length of alpha and omega. alpha
% is a vector of Fourier coefficients. omega is an array with two columns,
% where omega(k) is the frequnecy that corresponds to the coefficient
% alpha(k).
%
% The complexity of the computation is O(n*M^2).
%
% Yoel Shkolnisky, December 2009.

n=numel(alpha);

if size(omega,1)~=n
    error('alpha and omega must have the same number of rows');
end
if size(omega,2)~=2
    error('omega must have two columns')
end

low_idx=ceil((M-1)/2);
high_idx=floor((M-1)/2);
f=zeros(M);

for j1=-low_idx:high_idx
    for j2=-low_idx:high_idx
        for k=1:n
            v=alpha(k)*exp(2*pi*1i*(j1*omega(k,1)+j2*omega(k,2))/M);
            f(j1+low_idx+1,j2+low_idx+1)=f(j1+low_idx+1,j2+low_idx+1)+v;
        end
    end
end
