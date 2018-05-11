% VEC_TO_VOL Unroll vectors to volumes
%
% Usage
%    vol = vec_to_vol(vec);
%
% Input
%    vec: An N^3-by-... array.
%
% Output
%    vol: An N-by-N-by-N-by-... array.

function vol = vec_to_vol(vec)
    N = round(size(vec, 1)^(1/3));

    if N^3 ~= size(vec, 1)
        error('volumes must be cubic');
    end

    sz = size(vec);
    vol = reshape(vec, [N*ones(1, 3) sz(2:end)]);
end
