function r=Radial(w,org)
% function r = Radial(w,org)
% Compute the circularly-averaged, radial component of w,
% with the origin taken to be org.
% Given w is n x n, r is n/2 x 1.
% if w is three-dimensional, compute the averages over spherical shells.
% r(1) = average of points with r in [0,1),
% r(2) = average of points with r in [1,2), etc.
%
% Algorithm is based on that used in Hist.m

[nx ny nz] = size(w);
if nargin>1 || nx ~= ny || ndims(w)>2 % general version:
    if nargin<2
        org=[1 1 1]*ceil((nx+1)/2);
    end;
    if nz==1  % 2D input
        R=Radius([nx ny],org);
    else
        R=Radius3([nx ny nz],org);
    end;
    nr=floor(nx/2);
    r=zeros(nr,1);
    r0=zeros(nr,1);
    for i=1:nr
        q=R<i;
        r(i)=sum(q(:).*w(:));
        r0(i)=sum(q(:));
    end;
    r=[r(1);diff(r)];
%     r(2:nr)-r(1:nr-1);
    r0=[r0(1);diff(r0)];
%     =r0(2:nr)-r0(1:nr-1);
    r = r./r0;
else
    % fast Mex version
    [nx ny]=size(w);
    % we assume w is square and nx is even.
    R=Radius(nx);
    
    [r norm]=WeightedHisto(R,w,floor(nx/2));
    r=r./(norm+eps);
end;