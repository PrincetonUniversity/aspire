% VECMAT_TO_VOLMAT Roll up vector matrices into volume matrices
%
% Usage
%    volmat = vecmat_to_volmat(vecmat);
%
% Input
%    vecmat: A vector matrix of size N1^3-by-N2^3-by-... .
%
% Output
%    volmat: A volume "matrix" of size N1-by-N1-by-N1-by-N2-by-N2-by-
%       N2-by-... .

function volmat = vecmat_to_volmat(vecmat)
    sz = size(vecmat);

    N1 = round(sz(1)^(1/3));
    N2 = round(sz(2)^(1/3));

    if N1^3 ~= sz(1) || N2^3 ~= sz(2)
        error('volumes must be cubic');
    end

    volmat = reshape(vecmat, [N1*ones(1,3) N2*ones(1,3) sz(3:end)]);
end
