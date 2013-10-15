% function y=icfft(x)
%
% Aliased inverse Fourier transform (IFFT) of the sequence x.
% The FFT is computed using O(nlogn) operations.
%
% x    The sequence whose IFFT should be computed. 
%      Can be of odd or even length. Must be a 1D vector.
%
% Returns the aliased IFFT of the sequence x.
% 
% Yoel Shkolnisky 22/10/01

function y=icfft(x)
y = fftshift(ifft(ifftshift(x, 1), [], 1), 1);

% Revision Record
% 15/1/03	Yoel Shkolnisky		Use fftshift1d instead of fftshift