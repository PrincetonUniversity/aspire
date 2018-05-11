function d=SphereDensity(n,a,org)
% function d=SphereDensity(n,a,org)
% Compute the path length (i.e. integrated density) 
% through a sphere of radius a, projected onto an n x n image centered on
% the origin org.  By default org is ceil((n+1)/2)*[1 1]
% n x n image.  Each pixel is assigned a radius value r.  The value at that
% pixel is obtained as the integral of f(x)=sqrt(a^2-x^2) from r-1/2 to
% r+1/2.  In this way the resulting density d is smooth.

if numel(n)<2
    n=n*[1 1];
end;
if nargin<3
    org=ceil((n+1)/2);
end;
r=Radius(n,org);
rm=r-.5;
rp=r+.5;

d=real(rp.*sqrt(a^2-rp.^2)+a^2*atan(rp./(sqrt(a^2-rp.^2)))...
    -rm.*sqrt(a^2-rm.^2)-a^2*atan(rm./(sqrt(a^2-rm.^2))));

% imacs(d);
