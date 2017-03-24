function [ pf ]=cryo_pft_nfft(p, P)
%
% Compute the polar Fourier transform of projections with resolution n_r in
% the radial direction and resolution n_theta in the angular direction.
%
% If p is a volume, the function computes the polar Fourier transform of
% each slice in the volume seperately.
%
% Input parameters:
%   p          3D array.  p(:,:,k) is projection k.  
%              Two first dimensions are x and y of the each projection.
%   P          NFFT plan
%   P.freqs    frequency list for polar Fourier transform.
%   P.n_theta  Angular resolution. Number of Fourier rays computed for each
%              projection.
%   P.n_r      Number of samples on each radial line.
%
% Output parameters:
%   pf      Polar Fourier transform of the input array. pf is a matrix of
%           with n_r rows and n_theta columns. Each column corresponds to
%           a fixed theta. The first column corresponds a theta=0. The
%           last column corresponds to theta nearest 2*pi. The first row
%           corresponds to r=0. The lase row correspond to r nearest pi.
%           If f is a volume with n slices, pf is a volume of size 
%           n_r x n_theta x n. The third index is the slice number; the
%           other two are as above.
%       
% Zhizhen Zhao, 3/2015.
% Tejal Bhamre, 3/2017 : Using NFFT wrapper


% precomputed interpolation weights once for the give polar grid. This is
% used below for computing the polar Fourier transform of all slices


N1 = P.nL;
N2 = P.nL;
N=[N1;N2];
freqs = P.freqs;
M = length(freqs);

n_theta = P.n_theta;
n_r = P.n_r;
n_proj = size(p, 3);
pf=zeros(M, n_proj);

for i = 1:n_proj
    pf(:,i)=nufft2(p(:,:,i),-freqs'*2*pi);
end;

pf = reshape(pf, n_r, n_theta, n_proj);
