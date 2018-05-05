% BASIS_CHECK_EXPAND Checks input to basis expand functions
%
% Usage
%    basis_check_expand(x, basis);
%
% Input
%    x: The vector to be expanded in the basis.
%    basis: The basis object in which to expand.
%
% Description
%    Errors if the input does not correspond to the parameters of the basis
%    object.

function basis_check_expand(x, basis)
    sz_x = size(x);

    if numel(sz_x) < numel(basis.sz) || ...
        ~all(sz_x(1:numel(basis.sz)) == basis.sz)
        error('First dimensions of x must match basis.sz');
    end
end
