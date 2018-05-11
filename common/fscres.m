function r=fscres(fsc,cutoff)
%FSCRES  Find the resolution from an FSC curve.
%   r=fscres(fsc,cutoff)    Determine the index of the last bin of the
%   given fsc curve that is above the the given cutoff. Intermediate values
%   of the FSC curve are estimated using interpolation.
%
% Yoel Shkolnisky, December 2013.
%
% Revisions:
% Y.S March 2013    If the entire FSC curve is above the required cutoff,
%                   return the maxiaml possible resolution, that is, n.

n=numel(fsc);
r=n;
x=1:n;
xx=1:0.01:n;
y=interp1(x,fsc,xx,'PCHIP');
ii=find(y<cutoff,1,'first');
if ~isempty(ii)
    r=xx(ii);
end
