% FFT3 Calculate three-dimensional FFT
%
% Usage
%    x_f = fft3(x);
%
% Input
%    x: The three-dimensional signal to be transformed. The FFT is only
%       applied to the first three dimensions.
%
% Output
%    x_f: The standard Fourier transform of x.

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function x = fft3(x)
    x = fft(fft(fft(x, [], 1), [], 2), [], 3);
end
