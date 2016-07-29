% MAT_TO_VEC Converts a matrix into vectorized form
%
% Usage
%    vec = mat_to_vec(mat, is_symmat);
%
% Input
%    mat: An array of size N-by-N-by-... containing the matrices to be vector-
%       ized.
%    is_symmat: Specifies whether the matrices are symmetric/Hermitian, in
%       which case they are stored in packed form using symmat_to_vec (default
%       false).
%
% Output
%    vec: The vectorized form of the matrices, with dimension N^2-by-... or
%       N*(N+1)/2-by-... depending on the value of is_symmat.

function vec = mat_to_vec(mat, is_symmat)
    if nargin < 2 || isempty(is_symmat)
        is_symmat = false;
    end

    if ~is_symmat
        sz = size(mat);

        N = sz(1);

        if sz(2) ~= N
            error('matrix must be square');
        end

        vec = reshape(mat, [N^2 sz(3:end) 1]);
    else
        vec = symmat_to_vec(mat);
    end
end
