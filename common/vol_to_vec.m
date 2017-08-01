% VOL_TO_VEC Roll up volumes into vectors
%
% Usage
%    vec = vol_to_vec(vol);
%
% Input
%    vol: An N-by-N-by-N-by-... array.
%
% Output
%    vec: An N^3-by-... array.

function vec = vol_to_vec(vol)
    N = size(vol, 1);

    if size(vol, 2) ~= N || size(vol, 3) ~= N
        error('Volumes in `vol` must be cubic.');
    end

    sz = size(vol);
    vec = reshape(vol, [N^3 sz(4:end) 1]);
end
