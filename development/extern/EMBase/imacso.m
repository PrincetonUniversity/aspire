function imacso(x,y,m)
% function imacso(x,y,m);
% function imacso(m);
% Draw an image m whose center is at the origin; the
% same as imacs(fftshift(m)).
if nargin==3
    m=fftshift(m);
    imacs(x,y,m)
else
    imacs(fftshift(x));
end;
