function [noise,P,S,T1]=noise_exp2d(N,K,refpsd)
% NOISE_EXP2D   Generate colored noise with exponentially decaying auto
%               correlation.
% 
% [noise,P,S,T1]=noise_exp2d(N,K,noref)
% Generate K colored-noise images of size NxN.
% The isotropic autocorrelation of the noise is exp(-T1*r).
% T1 is set to 1 so that the autocovariance at distance 30 is about 9.e-14. 
%
% Input parameters:
%    N      Size of each noise image (NxN).
%    K      Number if noise images to generate.
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

% If one needs an autocorrelation that has the value delta at distance
% max_d, then the following code should be used:
% delta=1.0e-15;
% max_d=30;
% T1=-log(delta)/max_d; % Sampling rate

T1=1.66;
noise=zeros(N,N,K);
M=2*N-1;    % Need twice the number of random variables

[K1,K2]=meshgrid(-(N-1):(N-1),-(N-1):(N-1));
R=exp(-T1.*sqrt(K1.^2+K2.^2));
if norm(imag(R(:)))>1.0e-10
    warning('R has imaginary components');
end
H=real(cfft2(R)); % By Poisson formula, H is the periodized Fourier transform of exp(-r) times T1^2  
c2=M^2/(sum(H(:))); % The normalization ensures that the colored samples have variance 1
nf=ifftshift(sqrt(H.*c2));

W=iwindow(M,'hann');
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

% As the Fourier transform of R decays slowly, the filter computed above
% using cfft2(R) does not equal to the samples of the continuous Fourier
% transform. To compute the Fourier transform that corresponds to cfft2(R)
% we should us the Poisson summation formula.
%
% The samples of the true power spectrum are computed as follows:
% The function cryo_epsdR will estimate the continuous autocorrelation
% exp(-T1*r), whose Fourier transform 2*pi*T1/(T1^2+omega_r^2)^(3/2) will
% be sampled at the points 2*pi*k/(2*N+1).
%
S=-1;
if refpsd
    M=2*N-1; % T=M (T is sampling interval);
    omega0=2*pi/(2*N-1); % frequnecy steps
    omega1=omega0*M;
    omega_x=omega0.*K1;
    omega_y=omega0.*K2;
    c=2.*pi.*c2.*T1;
    S=zeros(size(P));    
    for j1=-100:100
        for j2=-100:100
%              omega_r=sqrt((omega_x+j1*omega1).^2+(omega_y+j2*omega1).^2);
%              S=S+2.*pi.*c2.*T1./((T1.^2+omega_r.^2).^(3/2)); % True power spectrum of the continuous process
            % Optimization:
            omega_r_sq=(omega_x+j1*omega1).^2+(omega_y+j2*omega1).^2;
            S=S+c./((T1.^2+omega_r_sq).^(3/2)); % True power spectrum of the continuous process
        end        
    end
    S=S./norm(S(:));
end
