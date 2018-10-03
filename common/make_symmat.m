% MAKE_SYMMAT Symmetrize a matrix
%
% Usage
%    B = make_symmat(A);
%
% Input
%    A: A matrix.
%
% Output
%    B: The Hermitian matrix (A+A')/2.

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function B = make_symmat(A)
    B = (A+A')/2;
end
