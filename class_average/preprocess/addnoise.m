function [noisy_projections, noise, I, sigma] = addnoise( projections, SNR, noise_type )
%Addnoise to projection images
%   noise_type: 'color' or 'gaussian'
%Zhizhen Zhao 09/01/2012

p = size(projections, 1);
K=size(projections, 3);
noisy_projections=zeros(size(projections));
noise=zeros(size(projections));
sigma=sqrt(var(reshape(projections(:, :, 1), p^2, 1))/SNR);
rand('state',1234);
randn('state',1234);
%for optimization, so we can use fft2 and
%ifft2 below instead of cfft2 and icfft2
if mod(p,2)==1
    lowidx=-(p-1)/2+p+1;
    highidx=(p-1)/2+p+1;
else 
    lowidx=-p/2+p+1;
    highidx=p/2+p;
end;

% Color Noise Response
I = cart2rad(2*p+1);
I=1./sqrt((1+I.^2));
%I = exp(-I/5);
I=I/norm(I(:));
noise_response = sqrt(I); 

for k=1:K
    gn=randn(2*p+1);
    if strcmpi(noise_type,'gaussian')
        cn=gn;
    else
        cn=real(icfft2(cfft2(gn).*noise_response));
    end
    cn=cn(lowidx:highidx,lowidx:highidx);
    cn=cn/std(cn(:));
    cn=cn.*sigma;
    noisy_projections(:,:,k)=projections(:, :, k) + cn;
    noise(:,:,k)=cn;
end

end

function [ I ] = cart2rad(N)
    N = floor(N);
    I = zeros(N,N);
  
    p = (N-1)/2;
    [X,Y] = meshgrid(-p:p,-p:p);
    
    I = sqrt(X.^2+Y.^2);
end

