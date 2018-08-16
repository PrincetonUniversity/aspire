% ISPREFIX Determines whether a vector is a prefix to another
%
% Usage
%    b = isprefix(x, y);
%
% Input
%    x, y: Two vectors.
%
% Output
%    b: True if `x` is shorter than `y` and `all(x == y(1:numel(x)))`, that is,
%       if `x` is a prefix of `y`.

function b = isprefix(x, y)
	if ~isvector(x)
		error('`x` must be a vector.');
	end

	if ~isvector(y)
		error('`y` must be a vector.');
	end

	b = (numel(x) <= numel(y) && all(x == y(1:numel(x))));
end
