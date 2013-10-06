function r=Radial2(w,org)
% function r = Radial2(w,org)
% Compute the circularly-averaged, radial component of the rectangular image w,
% with the origin taken to be org.  The number of points returned is
% nr=floor(min(size(w))/2).  For rectangular images r(nr) represents
% average values from the horizantal and vertical extreme points, so that
% Radial2(fftshift(abs(fftn(m)).^2)) gives the correct radial power
% spectrum for a rectangular image.
% r(1) = average of points with r in [0,1),
% r(2) = average of points with r in [1,2), etc.

sz = size(w);
if nargin<2
    org=ceil((sz+1)/2);
end;
szmin=min(sz);
nr=floor(szmin/2);
R=RadiusNorm(sz,org)*szmin;
    [r norm]=WeightedHisto(ceil(R),w,nr);
    r=r./(norm+eps);
return;
% 
% r=zeros(nr,1);
% r0=zeros(nr,1);
% for i=1:nr
%     p=find(R>=i-1 & R<i);
%     r(i)=sum(w(p));
%     r0(i)=numel(p);
% end;
% % r=[r(1);diff(r)];
% r = r./r0;

