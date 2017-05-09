function cryo_abinitio_C1_worker(algo,instack,outvol,outparams,showfigs,...
    verbose,n_theta,n_r,max_shift,shift_step)
% CRYO_ABINITO_C1_WORKER  Worker function for C1 abainitio reconstruction
%
% Internal function called by abinitio reconstruction functions.
% Do not call this function directly.
%
% Executes the C1 abinitio reconstruction algorithm specified by 'algo'.
% Takes as input an MRC stack of images and reconstructs an abinito model
% from the images. The reconstructed volume is stored in outvol.
% For a detailed description of the parameters see below.
%
% Parameters
%   algo        Orientation estimation algorithm to use. 1:sync3N,
%               2:sync2N, 3:LUD.
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

n_r_given=1;
if ~exist('n_r','var')
    n_r_given=0;
elseif n_r==-1
    n_r_given=0;
end

max_shift_given=1;
if ~exist('max_shift','var')
    max_shift_given=0;
elseif max_shift==-1
    max_shift_given=0;
end

if ~exist('shift_step','var')
    shift_step=1;
end

currentsilentmode=log_silent(verbose==0);

%% Load projections
projs=ReadMRC(instack);
K=size(projs,3);
log_message('projections loaded. Using K=%d projections of size %d x %d',K,size(projs,1),size(projs,2));
if size(projs,1)~=size(projs,2)
    error('Input images must be square');
end

%% Mask projections
mask_radius=round(size(projs,1)*0.45);
log_message('Masking projections. Masking radius is %d pixels',mask_radius);
[masked_projs,~]=mask_fuzzy(projs,mask_radius);

%% Compute polar Fourier transform
if ~n_r_given
    n_r=ceil(size(masked_projs,1)*0.5);
end
log_message('Computeing polar Fourier transform of projections. n_theta=%d, n_r=%d',n_theta,n_r);
[pf,~]=cryo_pft(masked_projs,n_r,n_theta,'single');  % take Fourier transform of projections

%% Find common lines from projections
if ~max_shift_given
    max_shift=ceil(size(projs,1)*0.15); % max_shift is 15% of the image size
end
log_message('Finding common lines using max_shift=%d, shift_step=%d',max_shift,shift_step);
[clstack,~,~,~]=cryo_clmatrix(pf,K,1,max_shift,shift_step);

log_message('Saving common lines');
save(outparams,'n_theta','n_r','clstack','max_shift','shift_step');

%% Orientation assigment
if algo==1
    % Use sync3N
    algname='sync3n';
    save(outparams,'algname','-append');
    
    VOTING_TICS_WIDTH=1;
    J_EIGS=4;
    J_WEIGHTS=true;
    S_WEIGHTS=true;
    
    % Estimate relative rotations
    [Rij0, r_valid_ks, r_good_ks, peak_width] = cryo_sync3n_estimate_all_Rijs...
        (clstack, n_theta, VOTING_TICS_WIDTH);    
    log_message('Stage 1 done: relative rotations');
    save(outparams,'Rij0','r_valid_ks','r_good_ks','peak_width','-append');
    
    % J-synchronization
    verbose = 2;
    [J_sync,J_significance,J_eigs,J_sync_iterations,~] = ...
        cryo_sync3n_Jsync_power_method (Rij0, J_EIGS, J_WEIGHTS, verbose);
    Rij = cryo_sync3n_flip_handedness(J_sync, Rij0);
    log_message('Stage 2 done: J-synchronization');
    save(outparams,'J_sync','J_significance','J_eigs','J_sync_iterations','-append');
    
    % Build 3NX3N Matrix S
    S = cryo_sync3n_syncmatrix(Rij);
    log_message('Stage 3 done: constructing matrix S');
    save(outparams,'S','-append');
    
    % S Weighting
    if S_WEIGHTS
        % Estimate weights for the 3x3 blocks of S
        [W, Pij, scores_hist] = cryo_sync3n_syncmatrix_weights(Rij0);
    else
        W = ones(N); % Default weights are all one (no weights)
        Pij = [];
        scores_hist = struct();
    end
    save(outparams,'S_WEIGHTS','W', 'Pij','scores_hist','-append');        
    log_message('Stage 4 done: computing weights');
    
    % Estimate rotations from S
    [rotations, S_eigs, ~] = cryo_sync3n_S_to_rot (S,10,W);
    save(outparams,'rotations','S_eigs','-append'); 
    log_message('Stage 5 done: estimating rotations');
    
    d_str=sprintf('%7.2f ',S_eigs);
    log_message('Top 10 eigenvalues of (weighted) sync matrix are %s',d_str);
    
elseif algo==2
    % Use sync2N
    algname='sync2n';
    save(outparams,'algname','-append');

    log_message('Starting buliding synchronization matrix');
    S=cryo_syncmatrix_vote(clstack,n_theta);
    log_message('Finished buliding synchronization matrix');
    save(outparams,'S','-append');
    
    rotations=cryo_syncrotations(S);
    save(outparams,'rotations','-append');
    
    d=eigs(S,10);
    d_str=sprintf('%7.2f ',d);
    log_message('Top 10 eigenvalues of sync matrix are %s',d_str);
    
elseif algo==3
    % Use LUD
    algname='LUD';
    save(outparams,'algname','-append');
    rotations = est_orientations_LUD(clstack, n_theta);
    save(outparams,'rotations','-append');
else
    error('Invalid vaule for ''algo''.');
end

if showfigs
    clerr=cryo_syncconsistency(rotations,clstack,n_theta);
    h=figure;
    hist(clerr(:,3),360);
end

log_message('Estimating shifts');
[est_shifts,~]=cryo_estimate_shifts(pf,rotations,max_shift,shift_step,10000,[],0);
save(outparams,'est_shifts','-append');
log_message('Finished estimating shifts');

%% Reconstruct downsampled volume with no CTF correction
n=size(projs,1);
[ v1, ~, ~ ,~, ~, ~] = recon3d_firm( projs,...
    rotations,-est_shifts, 1e-6, 100, zeros(n,n,n));
ii1=norm(imag(v1(:)))/norm(v1(:));
log_message('Relative norm of imaginary components = %e\n',ii1);
v1=real(v1);
WriteMRC(v1,1,outvol);

log_silent(currentsilentmode);
