% function y=cfft2(x)
%
% Aliased 2D FFT of the image x.
% The FFT is computed using O(n^2logn) operations.
%
% x   The image whose 2D FFT should be computed. Can be of odd or even length.
%
% Returns the aliased 2D FFT of the image x.
% 
% Yoel Shkolnisky 22/10/01


function y=cfft2(x)
	y = fftshift(fft2(ifftshift(x)));