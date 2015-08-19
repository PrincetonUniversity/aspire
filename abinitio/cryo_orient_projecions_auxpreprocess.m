function [vol,projs]=cryo_orient_projecions_auxpreprocess(vol,projs)
%
% Preprocess the volume and projections for the various
% cryo_orient_projections_* functions.
%
% Yoel Shkolnisky, August 2015.

szvol=size(vol);
% Downsample and normalize the given projection
if size(projs,3)==1
    projs=cryo_downsample(projs,[szvol(1) szvol(2)],0); % istack is zero
else
    projs=cryo_downsample(projs,[szvol(1) szvol(2)],1); % istack is zero
end
projs=cryo_normalize_background(projs,floor(szvol(1)/2),0);
