% IS_BASIS Checks whether object is a valid basis
%
% Usage
%    b = is_basis(basis);
%
% Input
%    basis: The object to be tested.
%
% Output
%    b: True is `basis` is a basis object, such as those returned by
%       `dirac_basis` or `fourier_bessel_basis`.

function b = is_basis(basis)
    b = true;

    b = b && isstruct(basis);

    required_fields = {'sz', 'count'};
    required_handles = {'evaluate', 'evaluate_t', 'expand', 'expand_t', ...
        'mat_evaluate', 'mat_evaluate_t', 'mat_expand', 'mat_expand_t'};

    b = b && all(isfield(basis, required_fields));
    b = b && all(isfield(basis, required_handles));

    b = b && ndims(basis.sz) == 2 && size(basis.sz, 1) == 1 && ...
        all(floor(basis.sz) == basis.sz) && all(basis.sz > 0);
    b = b && numel(basis.count) == 1 && ...
        floor(basis.count) == basis.count && basis.count > 0;

    check_handle = @(f)(isa(getfield(basis, f), 'function_handle'));
    b = b && all(cellfun(check_handle, required_handles));
end
