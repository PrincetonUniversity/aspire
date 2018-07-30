% MDIM_IFFTSHIFT Multi-dimensional FFT unshift
%
% Usage
%    x = mdim_ifftshift(x);
%    x = mdim_ifftshift(x, dims);
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
%    ifftshift, mdim_fftshift

function x = mdim_ifftshift(x, dims)
    if nargin < 2
        dims = 1:ndims(x);
    end

    for dim = dims
        x = ifftshift(x, dim);
    end
end
