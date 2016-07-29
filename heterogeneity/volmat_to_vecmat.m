% VOLMAT_TO_VECMAT Unroll volume matrices to vector matrices
%
% Usage
%    vecmat = volmat_to_vecmat(volmat);
%
% Input
%    volmat: A volume "matrix" of size N1-by-N1-by-N2-by-N2-by-N2-by-
%       N2-by-... .
%
% Output
%    vecmat: A vector matrix of size N1^3-by-N2^3-by-... .

function vecmat = volmat_to_vecmat(volmat)
    sz = size(volmat);

    if length(sz) < 6
        error('volume matrix must be at least six-dimensional');
    end

    N1 = sz(1);
    N2 = sz(4);

    if all(sz(1:6) ~= [N1*ones(1, 3) N2*ones(1, 3)])
        error('volumes must be cubic');
    end

    vecmat = reshape(volmat, [N1^3 N2^3 sz(7:end)]);
end
