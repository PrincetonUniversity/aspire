function g=nufft_t_2d_execute(beta,precomp)
%
% Compute the adjoint non-equally spaced FFT.
%
% The function computes the sums
%               n/2
%        g(j) = sum beta(k)*exp(i*k*x(j))
%              k=-n/2
% for j=1,...,n.
%
% The complexity of the algorithm is O(n*log(1/eps)+m*n(log(n)))
%
% The code is implemented such that eps is double precision. 
% m is the oversampling factor (m=2).
%
% Input parameters:
%    beta     Coefficients in the sums above. Real or complex numbers.
%    precomp  Precomputed interpolation factors, computed by the function
%             nufft_t_prepare_v3.
%
% This function is designed for computing the non-equally spaced of many
% images at the same frequnecies. Call the function nufft_t_2d_prepare
% once to compute the required data structure for the give frequnecies.
% 
% Note that x (the frequencies on which we sample the Fourier transform)
% should be longer than beta.
%
% Output:
%    g        The sums defined above.
%
% Yoel Shkolnisky, January 2008.
%
% Revisions:
% Filename changed from nufft_t_v3 to nufft_t_2d_prepare. Minor
% code fixes. Replace the call to nufftauxmx_v3 with a call to its renamed
% version nufftt2dexecutemx. (Y.S. December 22, 2009). 

if nargin<2
    error('Missing precomputed data');
end

b=precomp.b;
m=precomp.m;
q=precomp.q;

if mod(q,2)==1  % just for safety
    error('Set q to an even integer');
end

if length(size(beta))~=2
    error('beta must be a square matrix');
end

if size(beta,1)~=size(beta,2)
    error('beta must be a square matrix');
end

n=length(beta);
low_idx=-ceil((n-1)/2);
high_idx=floor((n-1)/2);
idx=low_idx:high_idx;

E=exp(b.*(2.*pi.*idx./(m*n)).^2);
E=E(:);
E=E*E.';
u=beta.*E;

low_idx_u=ceil((m*n-1)/2);
w=zeros(m*n,m*n);
w(low_idx_u+low_idx+1:low_idx_u+low_idx+1+n-1,low_idx_u+low_idx+1:low_idx_u+low_idx+1+n-1)=u;

W=fftshift(ifft2(ifftshift(w)));
W=W*numel(W);

g=nufftt2dexecutemx(W,precomp);