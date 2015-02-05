function [r theta]=RadiusNorm(n,org)
% function [r theta]=RadiusNorm(n,org)
% Create an n(1) x n(2) array where the value at (x,y) is the distance from the
% origin, normalized such that a distance equal to the width or height of
% the array = 1.  This is the appropriate function to define frequencies
% for the fft of a rectangular image.  
% For a square array of size n (or [n n]) the following is true:
% RadiusNorm(n) = Radius(n)/n.
% If a second output argument is given,
% the angle (in radians) is also supplied.
% The org argument is optional, in which case the FFT center is used.
% 
if numel(n)<2
    n=[n n];
end;
if nargin<2
     org=ceil((n+1)/2);
end;

% Commnted the following line. I want all numbers in double.
%n=single(n);  % this will force the output to be single
x=org(1);
y=org(2);
[X,Y]=ndgrid((1-x:n(1)-x)/n(1),(1-y:n(2)-y)/n(2));  % Make zero at x,y
    r=sqrt(X.^2+Y.^2);
if nargout>1
    theta=atan2(Y,X);
end;
