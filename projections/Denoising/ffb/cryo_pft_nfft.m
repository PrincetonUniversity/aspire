% CRYO_PFT_NFFT Compute polar Fourier transform
%
% Usage
%    pf = cryo_pft_nfft(proj, P);
%
% Input
%    proj: An array of images of size L-by-L-by-n.
%    P: The polar Fourier transform specification struct with the fields:
%           - freqs: The list of frequencies (see pft_freqs) for angles
%              between 0 and pi,
%           - n_theta: The number of samples in each concentric circle, and
%           - n_r: The number of samples in the radial direction.
%
% Output
%    pf: An array of size n_r-by-n_theta-by-n containing the non-uniform
%       Fourier transform at the frequncies specified by P.freqs.
%
% Description
%    Compute the polar Fourier transform of projections with resolution n_r in
%    the radial direction and resolution n_theta in the angular direction.
%
%    If p is a volume, the function computes the polar Fourier transform of
%    each slice in the volume separately.

% Written by Zhizhen Zhao - 3/2015.
% Modified by Tejal Bhamre (use NFFT wrapper) - 3/2017
% Reformatted and documented by Joakim Anden - 2018-Apr-13
% Optimized by Joakim Anden - 2018-Jul-24

function pf = cryo_pft_nfft(p, P)
    freqs = P.freqs;
    M = length(freqs);

    n_theta = P.n_theta;
    n_r = P.n_r;
    n_proj = size(p, 3);

    pf = nufft2(p, -freqs'*2*pi);

    pf = reshape(pf, n_r, n_theta/2, n_proj);
    pf = cat(2, pf, conj(pf));
end
