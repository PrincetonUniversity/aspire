function cryo_abinitio_C1_sync2n(instack,outvol,outparams,showfigs,...
    verbose,n_theta,n_r,max_shift,shift_step)
% CRYO_ABINITO_C1_SYNC2N  abinitio reconsturction using 3N-synchronization
%
% See cryo_abinitio_C1_sync3n for a detailed description of the paramters.
%
% Example:
% cryo_abinitio_C1_sync2n('projs.mrc','vol.mrc','molec.mat')
%
% Yoel Shkolnisky, March 2017.

% Check input and set default parameters
if ~exist('showfigs','var')
    showfigs=0;
end

if ~exist('verbose','var')
    verbose=1;
end

if ~exist('ntheta','var')
    n_theta=360;
end

if ~exist('n_r','var')
    n_r=-1;
end

if ~exist('max_shift','var')
    max_shift=-1;
end

if ~exist('shift_step','var')
    shift_step=1;
end

cryo_abinitio_C1_worker(2,instack,outvol,outparams,showfigs,...
    verbose,n_theta,n_r,max_shift,shift_step)