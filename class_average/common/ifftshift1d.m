% function y=ifftshift1d(x)
%
% Performs 1D ifftshift.
% A fast implementation of Matlab's ifftshift for 1D vectors.
%
% x   The vector to ifftshift. Must be a 1D vector.
%
% Yoel Shkolnisky 14/1/03

function y=ifftshift1d(x)

m=length(x);
p=floor(m/2);
y=x([p+1:m 1:p]);