% POSITIVE_FOURIER_MASK Get mask of positive Fourier frequencies
%
% Usage
%    [mask, real_mask] = positive_fourier_mask(sz);
%
% Input
%    sz: Size of the Fourier transform.
%
% Output
%    mask: A mask of those coefficients corresponding to "positive"
%       frequencies, in a lexicographical order.
%    real_mask: A mask of the coefficients that are real (that is, Nyquist-
%       type) within the positive mask.

function [mask, real_mask] = positive_fourier_mask(sz)
    for dim = 1:numel(sz)
        rngs{dim} = [0:floor(sz(dim)/2) -ceil(sz(dim)/2)+1:-1];
    end

    grids = cell(1, numel(sz));

    [grids{:}] = ndgrid(rngs{:});

    zero_mask = true(size(grids{1}));
    nonzero_mask = false(size(grids{1}));

    for dim = 1:numel(sz)
        nonzero_mask = nonzero_mask | ...
            (zero_mask&grids{dim}>0&grids{dim}<sz(dim)/2);
        zero_mask = zero_mask&(grids{dim}==0|grids{dim}==sz(dim)/2);
    end

    mask = zero_mask|nonzero_mask;
    real_mask = zero_mask;

    mask = mdim_fftshift(mask, 1:ndims(sz));
    real_mask = mdim_fftshift(real_mask, 1:ndims(sz));
end
