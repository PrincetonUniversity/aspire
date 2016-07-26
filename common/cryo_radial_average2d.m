function Pavg = cryo_radial_average2d(P)
% CRYO_RADIAL_AVERAGE2D Radial average of an image
%
% Pavg = cryo_radial_average2d(P)
%   Compute the radial average of the image P. The center is taken to be
%   the center of P. Pixel at distance d in Pavg is the average of all
%   pixels at distance d in P. P must be two-dimensional but need not be
%   square. The precision of Pavg is the same as the precision of P.
%
% Yoel Shkolnisky, August 2015.

if ~ismatrix(P)
    error('Input must be a 2D image');
end

sz=size(P);

if mod(sz(1),2)==1
    range1=-(sz(1)-1)/2:(sz(1)-1)/2;
else
    range1=-sz(1)/2+1/2:sz(1)/2-1/2;
end

if mod(sz(2),2)==1
    range2=-(sz(2)-1)/2:(sz(2)-1)/2;
else
    range2=-sz(2)/2+1/2:sz(2)/2-1/2;
end

[I,J]=ndgrid(range1,range2);

Pavg = zeros(size(P),class(P));
rsq=I.^2+J.^2;
rsqunique=unique(rsq);
[rsqsorted,rsqidx]=sort(rsq(:));
for i=1:numel(rsqunique)
    d=rsqunique(i);
     [li,ui]=bsearch(rsqsorted,d-1.0e-13,d+1.0e-13);
     Pavg(rsqidx(li:ui))=mean(P(rsqidx(li:ui)));
end;
