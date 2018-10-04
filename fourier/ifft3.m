% IFFT3 Calculate three-dimensional inverse FFT
%
% Usage
%    x = ifft3(x_f);
%
% Input
%    x_f: The three-dimensional signal to be transformed. The inverse FFT is
%       only applied to the first three dimensions.
%
% Output
%    x: The standard inverse Fourier transform of x_f.

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function x = ifft3(x)
    x = ifft(ifft(ifft(x, [], 1), [], 2), [], 3);
end
