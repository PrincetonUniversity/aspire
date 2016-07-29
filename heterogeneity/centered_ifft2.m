% CENTERED_FFT2 Calculate a centered, two-dimensional inverse FFT
%
% Usage
%    x = centered_ifft2(x_f);
%
% Input
%    x_f: The two-dimensional signal to be transformed. The inverse FFT is
%       only applied along the first two dimensions.
%
% Output
%    x: The centered inverse Fourier transform of x.

function x = centered_ifft2(x)
    x = ifftshift(ifftshift(x, 1), 2);
    x = ifft2(x);
    x = fftshift(fftshift(x, 1), 2);
end
