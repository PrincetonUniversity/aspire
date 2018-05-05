function [noise,P,S,T1]=noise_rexpr(N,K,refpsd)
% NOISE_REXPR   Generate colored noise with a decaying auto correlation.
% 
% [noise,P,S,T1]=noise_rexpr(N,K,noref)
% Generate K colored-noise images of size NxN.
% The isotropic autocorrelation of the noise is (1+r*T1)exp(-r*T1).
% For simplicity, T1=1.
%
% Input parameters:
%    N     Size of each noise image (NxN).
%    K     Number if noise images to generate.
%    refpsd Nonzero to compute reference psd. Default 0.
%
% Output parameters:
%    noise  K noise images of size NxN with the required autocorrelation.
%    P      Power spectrum of the noise estimated using the squared
%           absoulte value of the generated noise. 
%    S      Power spectrum of the noise. Array of size (2N-1)x(2N-1)
%           containing samples of the analytical Fourier transform at the
%           frequencies 2*pi/((2*N-1)*T1)*j, for j=-N:N.
%    T1     Sampling rate. Used for properly scaling the plots of the
%           estimated and analytic autocorrelations. 
%    
% Yoel Shkolnisky, May 2008.
%
% Revised: Yoel Shkolnisky, October 2014.

if ~exist('refpsd','var')
    refpsd=0;
end;

T1=1;
noise=zeros(N,N,K);
M=2*N-1;    % Need twice the number of random variables
[K1,K2]=meshgrid(-(N-1):N-1,-(N-1):N-1);
r=sqrt(K1.^2+K2.^2);
R=(1+r.*T1).*exp(-r.*T1);
if norm(imag(R(:)))>1.0e-10
    warning('R has imaginary components');
end
H=real(cfft2(R)); % By Poisson formula, H is the periodized Fourier transform of exp(-r) times T1^2  
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

% See noise_exp2d for a comment about computing the true power spectrum.
% Compute true power spectrum S using the Poisson formula - compute
% periodized power spectrum.
S=-1;
if refpsd   
    omega0=2*pi/(2*N-1); % frequnecy steps
    omega1=omega0*M; % bandlimit
    omega_x=omega0.*K1;
    omega_y=omega0.*K2;
    S=zeros(size(P));
    for k1=-100:100
        for k2=-100:100
            omega_x=omega_x+k1*omega1;
            omega_y=omega_y+k2*omega1;
            omega_r=sqrt(omega_x.^2+omega_y.^2);
            S=S+6*pi*T1^2./((T1^2+omega_r.^2).^(5/2));
        end
    end
    S=S./norm(S(:));
end

