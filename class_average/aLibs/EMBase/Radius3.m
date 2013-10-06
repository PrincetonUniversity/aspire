function r=Radius3(n,org,square)
% function r=Radius3(n,org,square)
% Create an n x n x n array where the value at (x,y,z) is the distance from the
% origin, the default origin being ceil(n/2+1,n/2+1,n/2+1).
% The org argument is optional, in which case the FFT center is used.
% If the square argument is given and is =1, r will be returned as the
% radius squared.
% The returned values are single precision.

if nargin<2
    org=floor([n/2+1 n/2+1 n/2+1]);
end;
if nargin<3
    square=0;
end;
% org=single(org);
n=single(n(1));  % this will force the output to be single
x=org(1);
y=org(2);
z=org(3);
[X,Y,Z]=ndgrid(1-x:n-x,1-y:n-y,1-z:n-z);  % Make zero at x,y,z
if square
    r=(X.^2+Y.^2+Z.^2);
else
    r=sqrt(X.^2+Y.^2+Z.^2);
end;
