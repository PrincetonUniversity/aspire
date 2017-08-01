% PRECOMPUTE_BASIS Precompute entire basis for fast evaluation/expansion
%
% Usage
%    basis = precompute_basis(basis);
%
% Input
%    basis: A basis object whose vectors are to be precomputed.
%
% Output
%    basis: The same basis object but with vectors precomputed.

function basis = precompute_basis(basis)
    B = basis.evaluate(eye(basis.count));

    B = reshape(B, [prod(basis.sz) basis.count]);

    basis.B = B;

    basis.evaluate = @(v)(precomputed_evaluate(v, basis));
    basis.expand = @(x)(precomputed_expand(x, basis));
    basis.evaluate_t = @(x)(precomputed_evaluate_t(x, basis));
    basis.expand_t = @(v)(precomputed_expand_t(v, basis));

    d = numel(basis.sz);

    basis.mat_evaluate = @(V)(mdim_mat_fun_conj(V, 1, d, basis.evaluate));
    basis.mat_expand = @(X)(mdim_mat_fun_conj(X, d, 1, basis.expand));
    basis.mat_evaluate_t = @(X)(mdim_mat_fun_conj(X, d, 1, basis.evaluate_t));
    basis.mat_expand_t = @(V)(mdim_mat_fun_conj(V, 1, d, basis.expand_t));
end

function x = precomputed_evaluate(v, basis)
    basis_check_evaluate(v, basis);

    [v, sz_roll] = unroll_dim(v, 2);

    x = basis.B*v;

    x = reshape(x, [basis.sz size(x, 2)]);

    x = roll_dim(x, sz_roll);
end

function v = precomputed_expand(x, basis)
    basis_check_expand(x, basis);

    [x, sz_roll] = unroll_dim(x, numel(basis.sz)+1);

    x = reshape(x, [prod(basis.sz) size(x, numel(basis.sz)+1)]);

    v = basis.B\x;

    v = roll_dim(v, sz_roll);
end

function v = precomputed_evaluate_t(x, basis)
    basis_check_expand(x, basis);

    [x, sz_roll] = unroll_dim(x, numel(basis.sz)+1);

    x = reshape(x, [prod(basis.sz) size(x, numel(basis.sz)+1)]);

    v = basis.B'*x;

    v = roll_dim(v, sz_roll);
end

function x = precomputed_expand_t(v, basis)
    basis_check_evaluate(v, basis);

    [v, sz_roll] = unroll_dim(v, 2);

    x = basis.B'\v;

    x = reshape(x, [basis.sz size(x, 2)]);

    x = roll_dim(x, sz_roll);
end
