% MDIM_MAT_FUN_CONJ Conjugate a multidimensional matrix using a linear mapping
%
% Usage
%    Y = mdim_mat_fun_conj(X, d1, d2, fun);
%
% Input
%    X: An N_1-by-...-by-N_d1-by-N_1...-by-N_d1-by-... array, with the first
%       2*d1 dimensions corresponding to matrices with columns and rows of
%       dimension d1.
%    d1, d2: The dimensions of the input (X) and output (Y) matrices,
%       respectively.
%    fun: A function handle of a linear map that takes an array of size
%        N_1-by-...-by-N_d1-by-... and returns an array of size
%        M_1-by-...-by-M_d2-by-... .
%
% Output
%    Y: An array of size M_1-by-...-by-M_d2-by-M_1-by-...-by-M_d2-by-...
%       resulting from applying fun to the rows and columns of the
%       multidimensional matrix X.

function Y = mdim_mat_fun_conj(X, d1, d2, fun)
    [X, sz_roll] = unroll_dim(X, 2*d1+1);

    X = fun(X);

    X = conj(permute(X, [d2+1:d2+d1 1:d2 d1+d2+1]));

    X = fun(X);

    X = conj(permute(X, [d2+1:2*d2 1:d2 2*d2+1]));

    X = roll_dim(X, sz_roll);

    Y = X;
end
