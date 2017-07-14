% FOURIER_FLIP Flip Fourier transform along specified dimensions
%
% Usage
%    yf = fourier_flip(xf, dims);
%
% Input
%    xf: The Fourier transform to be flipped.
%    dims: A vector of dimensions along which to flip.
%
% Output
%    yf: The Fourier transform xf, but flipped along the specified dimensions
%       with the zero frequency kept intact.

function xf = fourier_flip(xf, dims)
    idx_src.type = '()';
    idx_dst.type = '()';

    colon = {':'};

    idx_dst.subs = colon(ones(1, ndims(xf)));
    idx_src.subs = colon(ones(1, ndims(xf)));

    for dim = dims
        idx_dst.subs{dim} = 1:size(xf, dim);
        idx_src.subs{dim} = [1 size(xf, dim):-1:2];
    end

    xf = mdim_ifftshift(xf, dims);

    xf = subsasgn(xf, idx_dst, subsref(xf, idx_src));

    xf = mdim_fftshift(xf, dims);
end
