function out=cryo_mask(in,isstack,r,risetime)
% CRYO_MASK   Mask projections/volume
%
% out=cryo_mask(in)
%   Apply a mask to the volume 'in'. The volume must have equal dimensions.
%   If n is the side of the volume then the default radius of the mask is
%   floor(0.45*n), and the default risetime is floor(0.05*n).
%
% out=cryo_mask(in,isstack)
%   Set isstack to nonzero to interpreted the last dimension of 'in' as the
%   index of each image in the stack.
%
% out=cryo_mask(in,isstack,r)
%   Use masking radius 'r'.
%
% out=cryo_mask(in,isstack,r,risetime)
%   Use risetime 'risetime'.
%
% Yoel Shkolnisky, April 2017.

if ~exist('isstack','var')
    isstack=0;
end


nd=ndims(in);
if nd~=2 && nd~=3
    error('Input must be a 2D or 3D array');
end

if isstack
    nd=2;
end

sz=size(in);
if nd==2
    if sz(1)~=sz(2)
        error('Images of the input stack must have equal sides');
    end
elseif nd==3
    if ~all(sz==sz(1))
        error('Input volume must have equal sides');
    end
end

n=sz(1);
if ~exist('r','var')
    r=floor(0.45*n);
end

if ~exist('risetime','var')
    risetime=floor(0.05*n);
end
mask=fuzzymask(n,nd,r,risetime);
out=bsxfun(@times,in,mask);
