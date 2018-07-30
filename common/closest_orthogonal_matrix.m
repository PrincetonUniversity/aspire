% CLOSEST_ORTHOGONAL_MATRIX Project onto O(n)
%
% Usage
%    Q = closest_orthogonal_matrix(A);
%
% Input
%    A: An n-by-n matrix.
%
% Output
%    Q: An n-by-n orthogonal matrix that minimizes |A-Q|_F.

function Q = closest_orthogonal_matrix(A)
    [U, ~, V] = svd(A);

    Q = U*V';
end
