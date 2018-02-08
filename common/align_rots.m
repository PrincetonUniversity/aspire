% ALIGN_ROTS Align rotations matrices to a set of reference rotations
%
% Usage
%    rots_aligned = align_rots(rots, rots_ref);
%
% Input
%    rots: The rotations to be aligned in the form of a 3-by-3-by-n array.
%    rots_ref: The reference rotatons to which we would like to align in
%       the form of a 3-by-3-by-n array.
%
% Output
%    rots_aligned: The original rotations `rots`, optionally flipped and
%       rotated to align with `rots_ref`. Specifically, two sets of rotations
%       are computed, the original `rots` and their flipped copies, which
%       which correspond to conjugation by diag([1 1 -1]). These are then
%       rotated to minimize the Frobenius error with respect to `rots_ref` and
%       the closest set of rotations is returned.
%    Q: The 3-by-3 orthogonal matrix that aligns `rots` with `rots_ref`.
%    is_flipped: Specifies whether `rots` has to be flipped (conjugated by
%       diag([1 1 -1])) before alignment.
%
% See also
%    rot_error

function [rots_aligned, Q, is_flipped] = align_rots(rots, rots_ref)
    J = diag([1 1 -1]);

    rots_flip = mat_conj(rots, J);

    [err, Q] = rot_error(rots, rots_ref);
    [err_flip, Q_flip] = rot_error(rots_flip, rots_ref);

    if err_flip < err
        Q = Q_flip;
        rots = rots_flip;
        is_flipped = true;
    else
        is_flipped = false;
    end

    rots_aligned = submatfun(@(R)(Q*R), rots, 3);
end
