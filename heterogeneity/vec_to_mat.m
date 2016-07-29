% VEC_TO_MAT Converts a vectorized matrix into a matrix
%
% Usage
%    mat = vec_to_mat(vec, is_symmat);
%
% Input
%    vec: The vectorized representations. If the matrix is non-symmetric, this
%       array has the dimensions N^2-by-..., but if the matrix is symmetric,
%       the dimensions are N*(N+1)/2-by-... .
%    is_symmat: True if the vectors represent symmetric matrices. In this case
%       vec_to_symmat is called (default false).
%
% Output
%    mat: The array of size N-by-N-by-... representing the matrices.

function mat = vec_to_mat(vec, is_symmat)
    if nargin < 2 || isempty(is_symmat)
        is_symmat = false;
    end

    if ~is_symmat
        sz = size(vec);

        N = round(sqrt(sz(1)));

        if N^2 ~= sz(1)
            error('vector must represent square matrix');
        end

        mat = reshape(vec, [N*ones(1, 2) sz(2:end)]);
    else
        mat = vec_to_symmat(vec);
    end
end
