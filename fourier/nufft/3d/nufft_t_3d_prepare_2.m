function prepdata=nufft_t_3d_prepare_2(beta,precision)
%
% Compute nufft_t_3d for cases where we need to use the same volume many
% times with different samping points.
%
% In such cases, nufft_t_3d are replaced by 
%     prepdata=nufft_t_3d_prepare_2(beta,'single'); % or 'double'
%     g=nufft_t_3d_execute_2(x,prepdata);
% nufft_t_3d_execute_2 can be called many times with different x.
%
% See nufft_t_3d for information about the parameters.
%
% Yoel Shkolnisky, February 2010.

if nargin<2
    precision='single';
end

[b,m,q]=nufftparams(precision);

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
w = ifftshift(ifftshift(ifftshift(w, 1), 2), 3);
W = ifftn(w);
W = fftshift(fftshift(fftshift(W, 1), 2), 3);
W=W*numel(W);

prepdata.W=W;
prepdata.b=b;
prepdata.m=m;
prepdata.q=q;
prepdata.n=n;
