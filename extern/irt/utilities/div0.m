 function r = div0(num, den)
%function r = div0(num, den)
% r = num / den, except 0 where den = 0 
if nargin < 2, help(mfilename), error(mfilename), end

bad = den == 0;
r = num ./ (den + bad) .* (~bad);
