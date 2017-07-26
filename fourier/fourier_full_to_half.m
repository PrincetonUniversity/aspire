% FOURIER_FULL_TO_HALF Converts Fourier transform on full domain to half
%
% Usage
%    xf_half = fourier_full_to_half(xf, sig_sz);
%
% Input
%    xf: The Fourier transform defined on the full domain to be converted.
%    sig_sz: The size of the Fourier transform. This is equal to size(xf)
%       for dimensions 1 through d.
%
% Output
%    xf_half: The Fourier transform xf, restricted to its positive
%       frequencies as given by positive_fourier_mask. If size(xf, d+1) is
%       equal to n, xf will have n columns.

function xf_half = fourier_full_to_half(xf, sig_sz)
    d = numel(sig_sz);

    [xf, roll_sz] = unroll_dim(xf, d+1);

    n = size(xf, d+1);

    mask = positive_fourier_mask(sig_sz);

    idx = find(mask);

    xf = reshape(xf, [prod(sig_sz) n]);

    xf_half = xf(idx(:),:);

    xf_half = roll_dim(xf_half, roll_sz);
end
