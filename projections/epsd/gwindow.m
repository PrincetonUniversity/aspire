function w=gwindow(p,max_d)
% GWINDOW  Create 2D Gaussian window for spectral estimation.
%
% w=gwindpw(p)
% Return a (2p-1)x(2p-1) Gaussian window to be used for 2D power spectrum
% estimation. 
% 
% Input parameters:
%   p     Size of the returned mask.
%   max_d Width of the Gaussian.
%
% Output parameters:
%   w     (2p-1)x(2p-1) array.
%
% Yoel Shkolnisky, October 2014.

[X,Y] = meshgrid(-(p-1):(p-1),-(p-1):(p-1));
alpha=3.0; % Reciprocal of the standard deviation of the Gaussian window.
           % 1/alpha is the width of the Fourier transform of the window.
           % See Harris 78 for more details. 
w=exp(-alpha*(X.^2+Y.^2)/(2*max_d^2));
