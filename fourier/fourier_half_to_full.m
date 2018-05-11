% FOURIER_HALF_TO_FULL Converts Fourier transform on half domain to full
%
% Usage
%    xf = fourier_half_to_full(xf_half, sig_sz);
%
% Input
%    xf_half: An array of coefficicients corresponding to a Fourier transform
%       restricted to its positive frequencies as given by
%       positive_fourier_mask.
%    sig_sz: The full size of the Fourier transform.
%
% Output
%    xf: The full Fourier transform. If xf_half has n columns, xf will be of
%       size sig_sz-by-n.

function xf = fourier_half_to_full(xf_half, sig_sz)
    d = numel(sig_sz);

    [xf_half, roll_sz] = unroll_dim(xf_half, 2);

    n = size(xf_half, 2);

    [mask, real_mask] = positive_fourier_mask(sig_sz);

    idx = find(mask(:));
    real_idx = find(real_mask(:));

    xf = zeros([prod(sig_sz) n]);

    xf(idx(:),:) = xf_half;

    xf = reshape(xf, [sig_sz n]);

    xf_flip = fourier_flip(xf, 1:d);
    mask_flip = fourier_flip(mask, 1:d);

    xf = reshape(xf, [prod(sig_sz) n]);
    xf_flip = reshape(xf_flip, [prod(sig_sz) n]);

    xf(mask_flip,:) = conj(xf_flip(mask_flip(:),:));

    xf = reshape(xf, [sig_sz n]);

    xf = roll_dim(xf, roll_sz);
end
