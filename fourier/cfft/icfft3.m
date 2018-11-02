% ICFFT3 Calculate a centered, three-dimensional inverse FFT
%
% Usage
%    x = icfft3(x_f);
%
% Input
%    x_f: The three-dimensional signal to be transformed. The inverse FFT is
%       only applied along the first three dimensions.
%
% Output
%    x: The centered inverse Fourier transform of x_f.

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function x = icfft3(x)
    x = ifftshift(ifftshift(ifftshift(x, 1), 2), 3);
    x = ifft3(x);
    x = fftshift(fftshift(fftshift(x, 1), 2), 3);
end
