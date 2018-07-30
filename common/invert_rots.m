% INVERT_ROTS Invert an array of rotation matrices
%
% Usage
%    Qi = invert_rots(Q);
%
% Input
%    Q: A 3-by-3-by-n array of rotation matrices.
%
% Output
%    Qi: A 3-by-3-by-n array or rotation matrices such that each rotation
%       matrix is the inverse (that is, the transpose) of the input matrices Q.

function Qi = invert_rots(Q)
    Qi = permute(Q, [2 1 3]);
end
