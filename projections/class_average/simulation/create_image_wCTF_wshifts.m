function [proj, defocus_group, noise, noise_response, c, q]=create_image_wCTF_wshifts(data, SNR, noise_type)
%This function generates images with CTF with noise.
%Input: 
%       data.projections: simulated clean projection images. size: LxLxK
%       data.q: quaternions associated with each clean projection image.
%       size: 4xK
%       data.shifts: shifts for each image. matrix size :  Kx2
%       SNR: signal to noise ratio
%       noise_type: 'gaussian' additive white gaussian noise
%                   'color' colored noise with 1/sqrt(1+f^2) noise response
%                   'clean' clean images
%Output:
%       proj: projections images
%       defocus_group: indices for defocus groups
%       noise: noise images stack
%       sigma: variance of the noise
%       c: CTF functions
%       q: quaternions for rotation matrices
%
%Note this function needs precomputed dataset 70S_proj_10_5_s0
%Zhizhen Zhao Aug 2012


projections = data.projections;
q = data.q;
shifts = data.shifts;
clear data;
K = size(projections, 3);
L= size(projections, 1); %image size;
[projections]=shift_images(projections, shifts); %generate non-centered projection images.
%%CTF functions%%%%%%%%%%%%%%%%%%%%%%
lambda=EWavelength(200);
defocus=linspace(1.5, 4, 20);
num_defocus=length(defocus);
Cs=2.26;
B=100;
res=2.82;
alpha=0.07;
c=zeros(L, L, num_defocus);
for i=1:num_defocus
    c(:, :, i)=CTF(L, res, lambda, defocus(i), Cs, B, alpha);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize projections 
proj=zeros(size(projections));
a=floor(L/2);
defocus_group=zeros(K, 1);

for i=1:K
    tmp=cfft2(projections(:, :, i));
    tmp2=tmp.*c(:,:, mod(i, num_defocus)+1);
    proj(:, :, i)=real(icfft2(tmp2));
    defocus_group(i)=mod(i, num_defocus)+1;
end;

if ~strcmpi(noise_type,'clean')
    [proj, noise, noise_response]=addnoise(proj, SNR, noise_type);
end;

end
