% CFFT3 Calculate a centered, three-dimensional FFT
%
% Usage
%    x_f = cfft3(x);
%
% Input
%    x: The three-dimensional signal to be transformed. The FFT is only applied
%       along the first three dimensions.
%
% Output
%    x_f: The centered Fourier transform of x.

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function x = cfft3(x)
    x = ifftshift(ifftshift(ifftshift(x, 1), 2), 3);
    x = fft3(x);
    x = fftshift(fftshift(fftshift(x, 1), 2), 3);
end
