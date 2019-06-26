
% Find which pairs of rotations are in-plane rotations of each other
% Input: 
% rots = rotations
% filterAng= A pair of rotations which have beaming directions < filterAng
% are considered in-plane rotations of each other
% Output:
% inPlanePairs = nrot x nrot lower triangular matrix with non zero entries
% (ones) which mark pairs of rotations which are not in-plane rotations

function inPlanePairs=markInPlanePairs(rots,filterAng)
    projDirs=squeeze(rots(:,3,:));
    angles=projDirs'*projDirs;
    inPlanePairs=tril((abs(acos(angles))>filterAng)&...
        (abs(acos(angles))<(pi-filterAng)),-1);
end

