% function y=icfft2(x)
%
% Aliased 2D Inverse FFT of the sequence x.
% The inverse FFT is computed using O(n^2logn) operations.
%
% x    The frequency image whose 2D inverse FFT should be computed. Can be of odd or even length.
%
% Returns the aliased 2D inverse FFT of the sequence x.
% 
% Yoel Shkolnisky 06/02/02

function y=icfft2(x)
	y = fftshift(ifft2(ifftshift(x)));