% MAT_CONJ Conjugate matrix by another
%
% Usage
%    Y = mat_conj(X, M);
%
% Input
%    X: An N-by-N-by-... array where the first two indices correspond to
%       matrices.
%    M: An M-by-N matrix.
%
% Output
%    Y: An M-by-M-by-... array where the first two indices correspond to
%      the matrix products M*Z*M', where Z is a matrix formed by fixing all
%      indices three and higher in X.

function Y = mat_conj(X, M)
    Y = submatfun(@(Z)(M*Z*M'), X, 3);
end
