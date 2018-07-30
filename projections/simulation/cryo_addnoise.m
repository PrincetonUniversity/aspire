function [noisy_projections, noise, I, sigma] = cryo_addnoise( projections, SNR, noise_type,seed )
%
% Add additive noise to projection images.
%
% Input parameters:
%   projections     3D stack of projections. The slice projections(:,:,k)
%                   is the k'th projections.
%    SNR            Signal to noise of the output noisy projections.
%   noise_type      'color' or 'gaussian'
%   seed            Seed parameter for initializing the the random number
%                   generator of the noise samples. If not provided, the
%                   current state of the random number generator is used.
%                   In such a case the results of the function are not
%                   reproducible.
%
% Output parameters:
%    noisy_projections  Stack of noisy projections. Same size as the input
%                       projections.
%    noise              Stack containing the additive noise added to each
%                       projection. Same size as input projections.
%    I                  Normalized power spectrum of the noise. Has
%                       dimension equal to a single projection.
%    sigma              The standard deviation of Gaussian noise resulting
%                       in the required SNR. Computed using the first
%                       projection in the stack.
%
% Zhizhen Zhao 09/01/2012
% Revised:
% Y.S. August 2015  Add seed parameter.
%

p = size(projections, 1);
K=size(projections, 3);
noisy_projections=zeros(size(projections));
noise=zeros(size(projections));
sigma=sqrt(var(reshape(projections(:, :, 1), p^2, 1))/SNR);

%  Initialize the random number generator.
if exist('seed','var')
    initstate(seed); 
end

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
I1=ones(size(I));
I=1./sqrt((1+I.^2));
I=0*I1 + 1*I;
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
    p = (N-1)/2;
    [X,Y] = meshgrid(-p:p,-p:p);   
    I = sqrt(X.^2+Y.^2);
end

