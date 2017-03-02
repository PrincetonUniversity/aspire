function [resA,fighandle]=plotFSCA(vol1,vol2,cutoff,pixelsize,fighandle)
%PLOTFSC Draw Fourier shell correlation curve and estimate resolution.
%
% plotFSC(vol1,vol2,cutoff,pixelsize)
%   Align vol1 and vol2 and then draw the Fourier shell correlation curve
%   for the aligned volumes. Determine the resolution according the the
%   given cutoff and pixelsize. The volumes vol1 and vol2 need not be
%   aligned.
%
%   Input parameters:
%   vol1,vol2   3D arrays of the same dimensions, or names of MRC files.
%   cutoff  (Optional) Correlation cutoff threshold to determine
%           resolution. The resolution is determined by the first frequnecy
%           where the correlation goes below this value. Common values are
%           0.143 and 0.5. Default is 0.143.
%   pixelsize (Optional) Default=1A.
%   fighandle Where to plot the figure. If not given, a new handle is
%             allocated.
%
%   Output parameters:
%   resA    Resolution of the structure in Angstrom according to the given
%           cutoff value.
%   fighandle   Handle to the FSC curve plot.
%
% Example:
%   vol1=ReadMRC('./vol1.mrc');
%   vol2=ReadMRC('./vol2.mrc');
%   cutoff=0.143;
%   resA=plotFSCA(vol1,vol2,cutoff,1.5); 
%   fprintf('Reolsution from fsc50 curve is %5.2f\n',resA);
%
% See plotFSC for details on the algorithm for determining the resolution.
% Yoel Shkolnisky, May 2014.



if ~exist('pixelsize','var')
    pixelsize=1;
end

if ~exist('cutoff','var')
    cutoff=0.143;
end

if ~exist('fighandle','var')
    fighandle=figure;
end

%% Read volumes
if ~isnumeric(vol1)
    fname1=vol1;
    vol1=ReadMRC(fname1);
end

if ~isnumeric(vol2)
    fname2=vol2;
    vol2=ReadMRC(fname2);
end

%% Align vol1 and vol2
    [~,~,vol2aligned]=cryo_align_densities(vol1,vol2,pixelsize,1);

%% Plot FSC
[resA,fighandle]=plotFSC(vol1,vol2aligned,cutoff,pixelsize,fighandle);
