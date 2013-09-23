function precomp=nufft_t_3d_prepare_1(x,n,precision)
%
% Precomputataions for the 3D adjoint non-equally spaced FFT.
%
% Given a set of points x, the function generates all the required tables
% to compute the 3D adjoint non-equally spaced FFTs of volumes of size n^3.
% Call this function only once for given x, n, and precision, as
%       precomp=nufft_t_3d_prepare_1(x,n,precision)
% Then, for each volume beta, call
%       g=nufft_t_3d_execute_1(beta,precomp)
% beta must be of size nxnxn.
%
% Input parameters:
%    x          Sampling points. Real array with three columns of numbers
%               in the range [-pi,pi]^3. 
%    n          Size of the transformed volume.
%    precision  'double' or 'single'. Default is 'single'.
%    
%    x and the transformed volume (beta in nufft_t_3d_execute_1) need not
%    be the same length.
%
% Output:
%    precomp        Precomputed NUFFT tables.
%
% Yoel Shkolnisky, February 2010.

if nargin<3
    precision='single';
end

[b,m,q]=nufftparams(precision);

if length(size(x))~=2
    error('x must be a three-column array')
end

if size(x,2)~=3
    error('x must be a three-column array');
end

len_x=size(x,1);

nu=round(x*m*n/(2*pi));
beta1=zeros(len_x,q+1);
beta2=zeros(len_x,q+1);
beta3=zeros(len_x,q+1);

c1=(2*sqrt(b*pi))^3;
alpha = exp(-sum((x.*m.*n/(2*pi)-nu).^2,2)./(4*b))./c1;
gamma = exp(-((-q/2:q/2).').^2/(4*b));

for k=-q/2:q/2
    beta1(:,k+q/2+1)= exp(2.*(x(:,1).*m.*n/(2*pi)-nu(:,1)).*k/(4*b));
    beta2(:,k+q/2+1)= exp(2.*(x(:,2).*m.*n/(2*pi)-nu(:,2)).*k/(4*b));
    beta3(:,k+q/2+1)= exp(2.*(x(:,3).*m.*n/(2*pi)-nu(:,3)).*k/(4*b));
end

% Precompute the factors of Q
precomp.nu=nu;
precomp.alpha=alpha;
precomp.exp_k1_factor=beta1;
precomp.exp_k2_factor=beta2;
precomp.exp_k3_factor=beta3;
precomp.exp_k_square=gamma;
precomp.n=n;
precomp.m=m;
precomp.b=b;
precomp.q=q;



