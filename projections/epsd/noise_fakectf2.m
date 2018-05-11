function [noise,P,S]=noise_fakectf2(N,K)
% NOISE_FAKECTF2   Generate colored noise with known power spectrum.
% 
% [noise,P,S,T1]=noise_fakectf2(N,K)
% Generate K colored-noise images of size NxN.
% The isotropic power spectrum of the noise is supposed to simulate the
% true power spectrum of projections.
%
% Input parameters:
%    N     Size of each noise image (NxN).
%    K     Number if noise images to generate.
%
% Output parameters:
%    noise  K noise images of size NxN with the required autocorrelation.
%    P      Power spectrum of the noise estimated using the squared
%           absoulte value of the generated noise. 
%    S      Power spectrum of the noise. Array of size (2N-1)x(2N-1)
%           containing samples of the analytical Fourier transform.
%    
% Yoel Shkolnisky, May 2008.
%
% Revised: Yoel Shkolnisky, October 2014.

noise=zeros(N,N,K);
M=2*N-1;    % Need twice the number of random variables

sigma=1;
omega1=2*5*sigma;
omega0=omega1./M;

alpha=80*pi/omega1.^2;

% Sample the power spectrum 
[K1,K2]=meshgrid(-(N-1):N-1,-(N-1):N-1);
omega_x=omega0.*K1;
omega_y=omega0.*K2;
omega_r=sqrt(omega_x.^2+omega_y.^2);
H=(1./sqrt(2.*pi.*sigma.^2)).*exp(-omega_r.^2./(2*sigma.^2)).*(1+(sin(-alpha.*omega_r.^2)).^2);
c2=M^2/(sum(H(:))); % The normalization ensures that the colored samples have variance 1
nf=ifftshift(sqrt(H.*c2));

W=iwindow(M,'bartlett');
P=zeros(M,M,K);
for k=1:K % Apply filter
     gn=randn(M);              % Generate noise
     cn=ifft2(fft2(gn).*nf);   % Apply filter
     noise(:,:,k)=cn(1:N,1:N); % Take N samples
     P(:,:,k)=cfft2(cn.*W);
end

P=mean(abs(P).^2,3);
P=P./norm(P(:));


if norm(imag(noise(:)))>1.0e-8
    warning('Noise has imaginary components...');
end

noise=real(noise);

% Compute the reference power spectrum. 
omega0=omega1./(2*N-1);
omega_x=omega0.*K1;
omega_y=omega0.*K2;
omega_r=sqrt(omega_x.^2+omega_y.^2);
S=(c2./sqrt(2.*pi.*sigma.^2)).*exp(-omega_r.^2./(2*sigma.^2)).*(1+(sin(-alpha.*omega_r.^2)).^2);
S=S./norm(S(:));


