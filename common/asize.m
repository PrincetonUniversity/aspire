% ASIZE Size of an array
%
% Usage
%    sz = asize(x);
%
% Input
%    x: An array.
%
% Output
%    sz: The size of, excluding trailing 1s.
%
% Note
%    The behavior is the same as the builtin `size` function, except that it
%    returns a scalar for column vectors.

function sz = asize(x)
	sz = size(x);

	if numel(sz) == 2 && sz(2) == 1
		sz = sz(1);
	end
end
