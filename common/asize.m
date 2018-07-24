% ASIZE Size of an array
%
% Usage
%    sz = asize(x, dims);
%
% Input
%    x: An array.
%    dims: A set of dimensions to be retrieved (default `1:adims(x)`).
%
% Output
%    sz: The size of `x` along the dimensions specified by `dims`. If one of
%       `dims` is larger than `adims(x)`, that entry is assigned a value of
%       `1`.
%
% Description
%    Called with the default second argument `1:adims(x)`, the behavior is the
%    same as the builtin `size` function without second argument, except that
%    it returns a scalar for column vectors. When a second argument is given
%    it will behave the same as `size` when that argument is a scalar, with
%    the exception that `asize([], 2)` equals `1`, not `0`.

function sz = asize(x, dims)
	if nargin < 2
		dims = [];
	end

	n = max(1, adims(x));

	if isempty(dims)
		dims = 1:n;
	end

	sz0 = size(x);

	sz(dims > n) = 1;
	sz(dims <= n) = sz0(dims(dims <= n));
end
