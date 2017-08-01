% PAD_BASIS Zero-pad a basis to a certain size
%
% Usage
%    basis = pad_basis(basis, sz);
%
% Input
%    basis: The original basis object to be zero-padded.
%    sz: The desired size of the zero-padded basis. Note that the length of
%       sz must be the same as that of basis.sz to.
%
% Output
%    basis: A new basis object whose vectors are zero-padded to the desired
%       size.

function basis = pad_basis(basis, sz)
    if numel(basis.sz) ~= numel(sz)
        error('The number of elements in sz must be the same as basis.sz');
    end

    basis.type = [basis.type '-pad'];

    orig_sz = basis.sz;
    basis.sz = sz;

    basis.evaluate = @(v)(zero_pad(basis.evaluate(v), sz));
    basis.evaluate_t = @(x)(basis.evaluate_t(extract_center(x, orig_sz)));
    basis.expand = @(x)(basis.expand(extract_center(x, orig_sz)));
    basis.expand_t = @(v)(zero_pad(basis.expand_t(v), sz));

    d = numel(basis.sz);

    basis.mat_evaluate = @(V)(mdim_mat_fun_conj(V, 1, d, basis.evaluate));
    basis.mat_expand = @(X)(mdim_mat_fun_conj(X, d, 1, basis.expand));
    basis.mat_evaluate_t = @(X)(mdim_mat_fun_conj(X, d, 1, basis.evaluate_t));
    basis.mat_expand_t = @(V)(mdim_mat_fun_conj(V, 1, d, basis.expand_t));
end
