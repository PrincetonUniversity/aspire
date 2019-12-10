function r=Radial3(w,org)
% function r = Radial3(w,org)
% Compute the circularly-averaged, radial component of the cubic volume w,
% with the origin taken to be org.  The number of points returned is
% nr=floor(min(size(w))/2).
% r(1) = average of points with r in [0,1),
% r(2) = average of points with r in [1,2), etc.

sz = size(w);
if nargin<2
    org=[1 1 1]*ceil((sz+1)/2);
end;
szmin=min(sz);
nr=floor(szmin/2);
R=Radius3(sz,org);
    [r norm]=WeightedHisto(ceil(R),w,nr);
    r=r./(norm+eps);
return;
