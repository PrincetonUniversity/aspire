function masked_vol=cryo_mask_volume(vol,r,risetime)
%
% CRYO_MASK_VOLUME  Apply mask to a volume
%
% masked_vol=cryo_mask_volume(vol,r,risetime)
%   Mask the given volume with a circular mask made with an error
%   function with radius r and the given effective risetime (in pixels).
%   vol must be a 3D array with odd equal sides.
% masked_vol=cryo_mask_volume(vol,r)
%   Use risetime equal to 0.05 of the side of the volume.
% masked_vol=cryo_mask_volume(vol)
%   Use masking radius r which is 0.45 of the side of the volume, and
%   risetime which is 0.05 of the side of the volume.
%
% Yoel Shkolnisky, October 2016.

if ndims(vol)~=3
    error('Input must be a 3D array');
end

sz=size(vol);
if ~all(sz==sz(1))
    error('All dimensions of input array must be equal');
end

if mod(sz(1),2)==0
    error('Dimensions must be odd');
end

if ~exist('r','var')
    r=floor(0.45*sz(1));
end

if ~exist('risetime','var')
    risetime=floor(0.05*sz(1));
end

mask=fuzzymask(sz(1),3,r,risetime);
masked_vol=vol.*mask;