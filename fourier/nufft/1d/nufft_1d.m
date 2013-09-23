function f=nufft_1d(alpha,omega,precision,M)
%
% Non-equally spaced FFT.
%
% The function computes the sums
%                n
%        f(j) = sum alpha(k)*exp(2*pi*i*j*omega(k)/M)
%               k=1
% for j=-M/2,...,M/2, where omega(k) in [-n/2,n/2].
%
% The complexity of the algorithm is O(n*log(1/eps)+m*M(log(M)))
% m is the oversampling factor (m=2).
%
% The code was tested for M<n and everything seems to work fine.
%
% Input parameters:
%    alpha    Coefficients in the sums above. Real or complex numbers.
%    omega    Sampling frequnecies. Real numbers in the range [-n/2,n/2]. 
%    precision  'double' or 'single'. Default is 'single'.
%
% Output:
%    f        The sums defined above.
%
% Yoel Shkolnisky, December 2009.

if nargin<3
    precision='single';
end

[b,m,q]=nufftparams(precision);

n=numel(alpha);

if nargin<4
    M=n;
end

if length(size(alpha))~=2
    error('alpha must be a vector');
end

if (size(alpha,1)~=1) && (size(alpha,2)~=1)
    error('alpha must be a vector');
end

if length(size(omega))~=2
    error('omega must be a vector')
end

if size(omega,2)~=1
    error('omega must be a vector');
end

if size(omega,1)~=n
    error('alpha and omega must have the same number of rows');
end

mu=round(omega.*m);
P=zeros(n,q+1);
for j=-q/2:q/2
    P(:,j+q/2+1)=exp(-(m.*omega-(mu+j)).^2/(4*b))./(2*sqrt(b*pi));
end

tau=zeros(m*M,1);
offset1=ceil((m*M-1)/2);
for k=1:n
    for j=-q/2:q/2
        idx=mu(k)+j;
        idx=mod(idx+offset1,m*M);  
        tau(idx+1)=tau(idx+1)+P(k,j+q/2+1)*alpha(k);
    end
end

T=fftshift(ifft(ifftshift(tau)));
T=T*numel(T);

low_idx_M=-ceil((M-1)/2);
high_idx_M=floor((M-1)/2);
idx=low_idx_M:high_idx_M;
E=exp(b.*(2.*pi.*idx./(m*M)).^2);
E=E(:);
offset2=offset1+low_idx_M+1;
f=T(offset2:offset2+M-1).*E;

return;
