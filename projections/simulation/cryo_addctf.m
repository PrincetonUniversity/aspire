function ctfed_projections=cryo_addctf(projections,ctf)
% CRYO_ADDCTF Filter projections with a constrast transfer function (CTF).
%
%   ctfed_projections=CRYO_ADDCTF(projections,ctf) filters all projections
%   in the stack "projections" with the given CTF.
%
%   ctfed_projections=CRYO_ADDCTF(projections) uses a default CTF.
%
% ctf is constructed with the function CTF.
% Each projection in the stack "projections" must be suquare with odd
% sides.
% Output stack "ctfed_projections" has the same dimensions as the input
% stack.
%
% Example:
% K=5;
% dummySNR=1;
% projections= cryo_gen_projections(65,K,dummySNR);
% ctfed_projections=cryo_addctf(projections);
% figure;
% for k=1:K
%     % First row shows clean projections.
%     subplot(2,K,k); imagesc(projections(:,:,k)); 
%     axis image; axis off;
%     % Second row shows projections after CTF.
%     subplot(2,K,k+K); imagesc(ctfed_projections(:,:,k)); 
%     axis image; axis off;
% end
%
% Yoel Shkolnisky, September 2013.

nx=size(projections,1);
ny=size(projections,2);

if mod(nx,2)==0 || mod(ny,2)==0
    error('Only odd-sized images are supported');
end

if (nx~=ny)
    error('Projections must be square');
end
n=nx;

if ~exist('ctf','var')
    % Default CTF
    defocus=1.4;
    alpha=0.07;
    res=3.36;
    lambda=0.0251;
    Cs=2;B=100;
    ctf=CTF(n,res,lambda,defocus, Cs, B, alpha);
end  

if (size(ctf,1)~=n) || (size(ctf,2)~=n)
    error('CTF must have same dimensions as the proejctions');
end

K=size(projections,3);
ctfed_projections=zeros(size(projections));

for k=1:K
    p=projections(:,:,k);
    p_fourier=fft2(ifftshift(p));
    p_fourier=fftshift(p_fourier);
    p_fourier=p_fourier.*ctf; % Circular convolution is wrong. Fix this.
    p=fftshift(ifft2(ifftshift(p_fourier)));
    ctfed_projections(:,:,k)=p;
end
