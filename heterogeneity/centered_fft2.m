% CENTERED_FFT2 Calculate a centered, two-dimensional FFT
%
% Usage
%    x_f = centered_fft2(x);
%
% Input
%    x: The two-dimensional signal to be transformed. The FFT is only applied
%       along the first two dimensions.
%
% Output
%    x_f: The centered Fourier transform of x.

function x = centered_fft2(x)
    x = ifftshift(ifftshift(x, 1), 2);
    x = fft2(x);
    x = fftshift(fftshift(x, 1), 2);
end
