% ADIMS Number of dimensions excluding trailing singleton dimensions
%
% Usage
%    n = adims(x);
%
% Input
%    x: An array.
%
% Output
%    n: The number of dimensions of `x`, excluding trailing singleton
%       dimensions.
%
% Note
%    The behavior is the same as `ndims` except when `size(x, 2) == 1`, in
%    which case it returns `1` instead of `2`, or when `numel(x) == 0`, in
%    which case it returns `0`.

function n = adims(x)
	sz = size(x);

	n = numel(sz);

	if numel(sz) == 2 && sz(2) == 1 && sz(1) > 0
		n = 1;
	elseif numel(sz) == 2 && all(sz == 0)
		n = 0;
	end
end
