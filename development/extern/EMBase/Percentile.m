function [val upperVal]=Percentile(x,fraction)
%                 val=Percentile(x,fraction)
% [lowerVal upperVal]=Percentile(x,fraction)
% Returns the value of the element of x closest to the given fraction
% (between 0 and 1) of the distribution of x.  X can be of any size; all 
% elements are used.
% If two output variables are given, we find the lower and upper values
% corresponding to fraction and (1-fraction).  For example if fraction =
% .01, then the values of the 1% and 99% elements are returned.
x=x(:);
n=numel(x);
xs=sort(x);
nfrac=round(n*fraction)+1;
if nfrac<1 || nfrac>n+1
    error('Fraction out of bounds, must be in {0..1}');
end;
nfrac=min(n,nfrac);  % don't allow it to go beyond n
val=xs(nfrac);

if nargout>1
    upperNFrac=round(n*(1-fraction))+1;
    upperNFrac=min(n,upperNFrac);
    upperVal=xs(upperNFrac);
end;
    


