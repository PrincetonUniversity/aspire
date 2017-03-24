function Rij = cryo_sync3n_flip_handedness(J_sync, Rij)
% cryo_sync3n_flip_handedness   Flip handedness of relative rotations
%
% cryo_sync3n_flip_handedness(J_sync, Rij)
%   Take a list of relative rotations Rij and a vector J_sync of +/- 1 of
%   the same legnth, and flip the handedness of all relative rotations for
%   which J_sync is -1. J_sync can be computed using
%   cryo_sync3n_Jsync_power_method.
%   Returns the list of relative rotations Rij where rotations that
%   correspond to J_sync=-1 have been J conjugated.
%
% The function was revised from 
%   sync_J(J_sync, Rij)
%
% Ido Greenbeg, March 2016.
% Revised by Yoel Shkolnisky, March 2016.

J = diag([1 1 -1]);
for pair = 1:size(Rij,3)
    if J_sync(pair) < 0
        Rij(:,:,pair) = J*Rij(:,:,pair)*J;
    end
end

end
