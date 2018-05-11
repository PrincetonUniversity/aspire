% RAND_ROTS Generate random rotations
%
% Usage
%    rot_matrices = rand_rots(n);
%
% Input
%    n: The number of rotations to generate.
%
% Output
%    rot_matrices: An array of size 3-by-3-by-n containing n rotation matrices
%       sampled from the unifoorm distribution on SO(3).
%
% Note
%    This function depends on the random state of the `randn` function, so to
%    obtain consistent outputs, its state must be controlled prior to calling
%    using the `randn('state', s)` command.

function rot_matrices = rand_rots(n)
    qs = qrand(n);

    rot_matrices = q_to_rot(qs);
end
