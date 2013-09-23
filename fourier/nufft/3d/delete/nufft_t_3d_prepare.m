function prepdata=nufft_t_3d_prepare(beta,precision)
%
% Compute nufft_t_3d for cases where we need to use the same volume many
% times with different samping points.
%
% In such cases, nufft_t3_v2 and nufft_t3 are replaced by 
%     prepdata=nufft_t3_prepare(beta,'single'); % or 'double'
%     g=nufft_t3_execute(x,prepdata);
% nufft_t3_execute can be called many times with different x.
%
% See nufft_t3_v2 for information about the parameters.
%
% Yoel Shkolnisky, January 2009.
% Revisions:
%
% Filename changed from nufft_t3_prepare to nufft_t_3d_prepare. The call to
% nufft_t3_v2 was changed to a call to nufft_t_3d. (Y.S. December 22, 2009)
%

if nargin<2
    precision='single';
end

if (~strcmpi(precision,'single')) && (~strcmpi(precision,'double'))
    precision='single';
    warning('MATLAB:unkownOption','Unrecognized precsion. Using ''single''.');
end

if strcmpi(precision,'double')
    b=1.5629;
    m=2;
    q=28;
else %single precision
    b=0.5993;
    m=2;
    q=10;
end


if mod(q,2)==1
    error('Set q to an even integer');
end

if ndims(beta)~=3
    error('beta must be a cubed volume');
end

if (size(beta,1)~=size(beta,2)) || (size(beta,1)~=size(beta,2))
    error('beta must be a cubed volume');
end

n=size(beta,1);

low_idx=-ceil((n-1)/2);
high_idx=floor((n-1)/2);
idx=low_idx:high_idx;

E=exp(b.*(2.*pi.*idx./(m*n)).^2);
E=E(:);

EE=E*E.';
E3=zeros(n,n,n);
for k=1:n
    E3(:,:,k)=EE.*E(k);
end

u=beta.*E3;
low_idx_u=ceil((m*n-1)/2);
w=zeros(m*n,m*n,m*n);
w(low_idx_u+low_idx+1:low_idx_u+low_idx+1+n-1,...
  low_idx_u+low_idx+1:low_idx_u+low_idx+1+n-1,...
  low_idx_u+low_idx+1:low_idx_u+low_idx+1+n-1)=u;

% Compute the 3D FFT of  the padded volume only once, and store it in
% prepdata.
W=fftshift(ifftn(ifftshift(w)));
W=W*numel(W);

prepdata.W=W;
prepdata.b=b;
prepdata.m=m;
prepdata.q=q;
prepdata.n=n;