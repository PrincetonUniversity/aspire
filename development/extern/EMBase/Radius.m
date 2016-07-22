function [r theta]=Radius(n,org,square)
% function [r theta]=Radius(n,org,square)
% Create an n(1) x n(2) array where the value at (x,y) is the distance from the
% origin, the default origin being ceil((n+1)/2).  If n is a scalar,
% then the output array is n x n.
% If a second output argument is given,
% the angle (in radians) is also supplied.
% The org argument is optional, in which case the FFT center is used.
% If the square argument is given and is =1, r will be returned as the
% radius squared.
% The returned values are single precision.

if numel(n)<2
    n=[n n];
end;
if nargin<2
     org=ceil((n+1)/2);
end;
if nargin<3
    square=0;
end;

n=single(n);  % this will force the output to be single
x=org(1);
y=org(2);
[X,Y]=ndgrid(1-x:n(1)-x,1-y:n(2)-y);  % Make zero at x,y
if square
    r=(X.^2+Y.^2);
else
    r=sqrt(X.^2+Y.^2);
end;

if nargout>1
    theta=atan2(Y,X);
end;
