% DIRAC_BASIS Construct a Dirac basis object
%
% Usage
%    basis = dirac_basis(sz, mask);
%
% Input
%    sz: The size of the vectors for which to define the basis.
%    mask: A boolean mask of size sz indicating which coordinates to include
%       in the basis (default true(sz)).
%
% Output
%    basis: A Dirac basis object corresponding to the parameters.
%
% Description
%    The returned basis object can be used to project a vector x in the
%    standard basis onto the basis or to map a coefficient vector v into
%    the standard basis. This is achieved using the basis.expand and
%    basis.evaluate functions, respectively. The fields basis.sz denotes
%    the size of the vectors in the standard basis while basis.count denotes
%    the number of vectors in the basis.
%
%    For example, we can generate a random vector in the standard basis and
%    project it onto the Dirac basis through
%
%       basis = dirac_basis(sz, mask);
%       x = randn([sz 1]);
%       v = basis.expand(x);
%
%    Likewise, we can map a random coefficient vector into the standard basis
%    using
%
%       v = randn(basis.count, 1);
%       x = basis.evaluate(v);

function basis = dirac_basis(sz, mask)
    if nargin < 2 || isempty(mask)
        mask = true(sz);
    end

    basis = struct();

    basis.type = 'dirac';

    basis.sz = sz;
    basis.count = sum(mask(:));

    basis.mask = mask;

    basis.evaluate = @(v)(dirac_evaluate(v, basis));
    basis.expand = @(x)(dirac_expand(x, basis));

    basis.evaluate_t = basis.expand;
    basis.expand_t = basis.evaluate;

    d = numel(basis.sz);

    basis.mat_evaluate = @(V)(mdim_mat_fun_conj(V, 1, d, basis.evaluate));
    basis.mat_expand = @(X)(mdim_mat_fun_conj(X, d, 1, basis.expand));

    basis.mat_evaluate_t = @(X)(mdim_mat_fun_conj(X, d, 1, basis.evaluate_t));
    basis.mat_expand_t = @(V)(mdim_mat_fun_conj(V, 1, d, basis.expand_t));
end

function x = dirac_evaluate(v, basis)
    basis_check_evaluate(v, basis);

    [v, sz_roll] = unroll_dim(v, 2);

    x = zeros([prod(basis.sz) size(v, 2)], class(v));

    x(basis.mask,:) = v;

    x = reshape(x, [basis.sz size(x, 2)]);

    x = roll_dim(x, sz_roll);
end

function v = dirac_expand(x, basis)
    basis_check_expand(x, basis);

    [x, sz_roll] = unroll_dim(x, numel(basis.sz)+1);

    x = reshape(x, [prod(basis.sz) size(x, numel(basis.sz)+1)]);

    v = x(basis.mask,:);

    v = roll_dim(v, sz_roll);
end
