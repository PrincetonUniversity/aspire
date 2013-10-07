% function y=icfftn(x)
%
% Aliased n-dimensional Inverse FFT of the array x.
% The inverse FFT is computed using O((n^d)logn) operations, 
% where d is the dimension of the image.
%
% x    The frequency image whose inverse FFT should be computed.
%      Can be of odd or even length in each dimension.
%
% Returns the aliased n-dimensional inverse FFT of the array x.
% 
% Yoel Shkolnisky 25/02/03

function y=icfftn(x)
	y = fftshift(ifftn(ifftshift(x)));