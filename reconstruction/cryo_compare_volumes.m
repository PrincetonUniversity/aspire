function [resA,vol_aligned,h]=cryo_compare_volumes(vol,emd_or_vol,cutoff,pixA,verbose)
%
% CRYO_COMPARE_VOMLUES  Estimate resolution of volumes
% 
% [fsc,vol_aligned,h]=cryo_compare_volumes(vol,emd_or_vol)
%   Compute the the resolution of two volumes using the given cutoff
%   criterion of the Fourier shell corelation.
%   The second volume (emd_or_vol) is downsampled to have the same
%   dimensions as vol.
%
% Input parameters
%   vol  3D array or a filename in MRC format.
%   emd_or_vol Numeric value with the EMDID of the reference volume, or a
%              filename, or a 3D array.
%   cutoff  Cutoff value of the Fourier shell correlation to determine the
%           resolution. Default: 0.143.
%   pixA    Pixel size in Angstrom.
%   verbose Set to nonzero to print verbose messages.
%
% Output parameters:
%   resA    Agreement resolution of the two volumes.
%   vol_aligned     vol after alignment to emd_or_vol.
%   h   handle of the Fourier shell correlation plot.
%
% Yoel Shkolnisky, January 2018.

if ~exist('verbose','var')
    verbose=0;
end

if ~exist('pixA','var')
    pixA=1;
end

if ~exist('cutoff','var')
    cutoff=0.143;
end

if isscalar(emd_or_vol) % EMDID was provided
    mapfile=cryo_fetch_emdID(emd_or_vol,verbose);
    volref=ReadMRC(mapfile);
elseif ischar(emd_or_vol)
    volref=ReadMRC(emd_or_vol);
else
    volref=emd_or_vol;   
end

if ischar(vol)
    vol=ReadMRC(vol);
end

% Downsample volref to the dimensions of vol
sz=size(vol);
if verbose
    log_message('Downsampling volref to size %dx%dx%d',sz(1),sz(2),sz(3));
end
volref=cryo_downsample(volref,sz,0);
% Align vol to volref
if verbose
    log_message('Aligning volumes');
end
%[~,~,vol_aligned]=cryo_align_densities(volref,vol,pixA,verbose,cutoff);
[~,~,~,vol_aligned]=cryo_align_vols(volref,vol,verbose);

% Plot FSC
n=size(volref,1);
volref_masked=volref.*fuzzymask(n,3,floor(0.45*n),floor(0.05*n));
vol_aligned_masked=vol_aligned.*fuzzymask(n,3,floor(0.45*n),floor(0.05*n));

[resA,h]=plotFSC(volref_masked,vol_aligned_masked,cutoff,pixA);
if verbose
    log_message('Resolution between masked volumes is %5.2fA (using cutoff %4.3f)',resA,cutoff);
end
