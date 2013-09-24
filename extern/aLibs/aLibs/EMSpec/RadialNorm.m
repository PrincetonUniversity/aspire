function r=RadialNorm(w,nbins,org)
% function r = RadialNorm(w,nbins,org)
% Compute the circularly-averaged, radial component of w,
% with the origin taken to be org, and the maximum radius corresponding to
% nbins.  Thus if size(w)=[nx ny] then the farthest bin will correspond to
% the values in an ellipse with radii nx/2 and ny/2.
% For example, to get the circularly-averaged spectrum of an image use
%  S=RadialNorm(fftshift(abs(fftn(m)).^2),nbins);
% Only works in 2d at the moment.
n=size(w);
if numel(n)<2
    n(2)=n;
end;
if nargin<2
    nbins=floor(min(n/2));
end;
if nargin<3
    org=ceil((n+1)/2);
end;
% we assume w is square and nx is even.
R=RadiusNorm(n,org);

[r norm]=WeightedHisto(R*nbins*2,w,nbins);
r=r./(norm+eps);
