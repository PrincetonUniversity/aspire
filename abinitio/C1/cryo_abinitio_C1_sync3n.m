function cryo_abinitio_C1_sync3n(instack,outvol,outparams,showfigs,...
    verbose,n_theta,n_r,max_shift,shift_step,map_filter_radius)
% CRYO_ABINITO_C1_SYNC3N  abinitio reconsturction using 3N-synchronization
%
% Parameters
%   instack     Name of MRC file containing the projections (or class
%               averages) from which to estimate an abinitio model.
%   outvol      Name of MRC file into which to save the reconstructed
%               volume.
%   outparams   Name of MAT file in which intermediate outcomes of the
%               reconstruction algorithm are save (such as estimated
%               rotations and shifts). Used to debugging and detailed
%               analysis of the results.
%   showfigs    (optional) Set to 1 to display quaility assessment figures.
%               Currently only consistency histogram for the estimated 
%               rotations is display. Default is 0.
%   verbose     (Optional) Set to 1 to print verbose long messages. Default
%               is 1.
%   ntheta      (Optional) Angular resolution for common lines detection.
%               Default 360. 
%   n_r         (Optional) Radial resolution for common line detection.
%               Default is half the width of the images.
%   max_shift   (Optional) Maximal 1d shift (in pixels) to search between
%               common-lines. Default is 15% of image width of the images.
%   shift_step  (Optional) Resolution of shift estimation in pixels. Note
%               that shift_step can be any positive real number. Default:1. 
%   map_filter_radius   Detect common lines bet similar nighborhoods in the
%               images and not just similar Fourier rays. See
%               cryo_clmatrix_cpu for more details.
%
% Setting a parameter to empty value results in using the default value for
% this parameter.
%
% Example:
% cryo_abinitio_C1_sync3n('projs.mrcs','vol.mrc','molec.mat')
%
% Yoel Shkolnisky, March 2017.

% Check input and set default parameters
if ~exist('showfigs','var') || isempty(showfigs)
    showfigs=0;
end

if ~exist('verbose','var') || isempty(verbose)
    verbose=1;
end

if ~exist('n_theta','var') || isempty(n_theta)
    n_theta=360;
end

if ~exist('n_r','var') || isempty (n_r)
    n_r=-1;
end

if ~exist('max_shift','var') || isempty(max_shift)
    max_shift=-1;
end

if ~exist('shift_step','var') || isempty(shift_step)
    shift_step=1;
end

if ~exist('map_filter_radius','var') || isempty(map_filter_radius)
    map_filter_radius=0;
end

cryo_abinitio_C1_worker(1,instack,outvol,outparams,showfigs,...
    verbose,n_theta,n_r,max_shift,shift_step,map_filter_radius)
