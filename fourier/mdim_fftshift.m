% MDIM_FFTSHIFT Multi-dimensional FFT shift
%
% Usage
%    x = mdim_fftshift(x);
%    x = mdim_fftshift(x, dims);
%
% Input
%    x: The array to be shifted.
%    dims: An array of dimension indices along which the shift should occur.
%       If left empty, the shift is performed along all dimensions.
%
% Output
%    x: The x array shifted along the desired dimensions.
%
% See also
%    fftshift, mdim_ifftshift

function x = mdim_fftshift(x, dims)
    if nargin < 2
        dims = 1:ndims(x);
    end

    for dim = dims
        x = fftshift(x, dim);
    end
end
