% function y=fftshift1d(x)
%
% Performs 1D fftshift.
% A fast implementation of Matlab's fftshift for 1D vectors.
%
% x   The vector to fftshift. Must be a 1D vector.
%
% Yoel Shkolnisky 14/1/03

function y=fftshift1d(x)

m=length(x);
p=ceil(m/2);
y=x([p+1:m 1:p]);