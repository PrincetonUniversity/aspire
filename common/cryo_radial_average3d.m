function Vavg = cryo_radial_average3d(V)
% CRYO_RADIAL_AVERAGE3D Radial average of a volume
%
% Pavg = cryo_radial_average3d(V)
%   Compute the radial average of the volume P. The center is taken to be
%   the center of P. Pixel at distance d in Pavg is the average of all
%   pixels at distance d in P. P must be two-dimensional but need not be
%   square. The precision of Pavg is the same as the precision of P.
%
% Yoel Shkolnisky, August 2015.

if ndims(V)~=3
    error('Input must be a 3D volume');
end

sz=size(V);

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

if mod(sz(3),2)==1
    range3=-(sz(3)-1)/2:(sz(3)-1)/2;
else
    range3=-sz(3)/2+1/2:sz(3)/2-1/2;
end


[X,Y,Z]=ndgrid(range1,range2,range3);

Vavg = zeros(size(V),class(V));
rsq=X.^2+Y.^2+Z.^2;
rsqunique=unique(rsq);
[rsqsorted,rsqidx]=sort(rsq(:));
for i=1:numel(rsqunique)
    d=rsqunique(i);
     [li,ui]=bsearch(rsqsorted,d-1.0e-13,d+1.0e-13);
     Vavg(rsqidx(li:ui))=mean(V(rsqidx(li:ui)));
end;
