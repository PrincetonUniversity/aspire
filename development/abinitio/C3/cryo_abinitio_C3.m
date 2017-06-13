function cryo_abinitio_C3(instack,outvol,outparams,...
    n_projs,n_theta,n_r,max_shift,shift_step)
% CRYO_ABINITO_C4  abinitio reconsturction of a c4 symmetric molecule
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
%   ntheta      (Optional) Angular resolution for common lines detection.
%               Default 360. 
%   n_projs     (Optional) Number of projections to use from 'instack'. -1
%               means to use all projections. Default: -1 (all projections)
%   n_r         (Optional) Radial resolution for common line detection.
%               Default is half the width of the images.
%   max_shift   (Optional) Maximal 1d shift (in pixels) to search between
%               common-lines. Default is 15% of image width of the images.
%   shift_step  (Optional) Resolution of shift estimation in pixels. Note
%               that shift_step can be any positive real number. Default:1. 
%
% Example:
% cryo_abinitio_C4('/mnt/ix2/backup/datasets/10005/output/averages_nn100_group1.mrc','vol.mrc','molec_c4.mat')
%

% Check input and set default parameters
if ~exist('n_theta','var')
    n_theta = 360;
end

if ~exist('n_projs','var')
    n_projs_given = false;
else
    n_projs_given = true;
end

if ~exist('n_r','var')
    % defer exact value once the size of images is known
    n_r_given = false;
else
    n_r_given = true;
end

if ~exist('max_shift','var')
    % defer maximum shift once the size of images is known
    max_shift_given = false;
else
    max_shift_given = true;
end

if ~exist('shift_step','var')
    shift_step = 0.5;
end

%% Load projections
projs = ReadMRC(instack);

if n_projs_given
    if n_projs == -1
        nImages = size(projs,3);
    elseif n_projs <= 0 
        error('n_projs must be either positive number, or -1 for using all images');
    else
        assert(n_projs <= size(projs,3));
        nImages = n_projs;
    end
else % not provided as a parameter so use everything
    nImages = size(projs,3);
end

% im_indeces = floor(linspace(1,size(projs,3),nImages));
im_indeces = randperm(size(projs,3),nImages);
save(outparams,'im_indeces');

projs = projs(:,:,im_indeces);
assert(size(projs,3) == nImages);

log_message('projections loaded. Using %d projections of size %d x %d',nImages,size(projs,1),size(projs,2));
if size(projs,1)~=size(projs,2)
    error('Input images must be square');
end

%% Mask projections
mask_radius = round(size(projs,1)*0.4);
log_message('Masking projections. Masking radius is %d pixels',mask_radius);
[masked_projs,~] = mask_fuzzy(projs,mask_radius);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 1  : Computing polar Fourier transform of projections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~n_r_given
    n_r = ceil(size(masked_projs,1)*0.5);
end

[npf,~] = cryo_pft(masked_projs,n_r,n_theta,'single');

npf = gaussian_filter_imgs(npf);

save(outparams,'n_theta','n_r');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 2  : detect a single pair of common-lines between each pair of images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~max_shift_given
    max_shift = ceil(size(projs,1)*0.15); % max_shift is 15% of the image size
end
log_message('detecting common-lines');
clmatrix = cryo_clmatrix(npf,nImages,1,max_shift,shift_step); 
save(outparams,'clmatrix','max_shift','shift_step','-append');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 3  : detect self-common-lines in each image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sclmatrix = cryo_self_clmatrix_gpu(npf,max_shift,shift_step);
save(outparams,'sclmatrix','-append');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 4  : calculate self-relative-rotations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Riis = estimate_all_Riis(sclmatrix,n_theta);
save(outparams,'Riis','-append');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 5  : calculate relative-rotations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rijs = cryo_c3_estimate_all_Rijs(clmatrix,n_theta);
save(outparams,'Rijs','-append');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 6  : inner J-synchronization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[vijs,viis,im_inds_to_remove,pairwise_inds_to_remove,...
    npf,projs] = local_sync_J(Rijs,Riis,npf,projs);
save(outparams,'vijs','viis','im_inds_to_remove','pairwise_inds_to_remove','-append');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 7  : outer J-synchronization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[vijs,viis] = global_sync_J(vijs,viis);
save(outparams,'vijs','viis','-append');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 8  : third rows estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vis  = estimate_third_rows(vijs,viis);
save(outparams,'vis','-append');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 9  : in-plane rotations angles estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rots = estimate_inplane_rotations2(npf,vis,1,max_shift,shift_step);
save(outparams,'rots','-append');

estimatedVol = reconstruct(projs,rots,n_r,n_theta);   
WriteMRC(estimatedVol,1,outvol);
% 
% 
% log_message('Estimating shifts');
% [est_shifts,~]=cryo_estimate_shifts(pf,rotations,max_shift,shift_step,10000,[],0);
% save(outparams,'est_shifts','-append');
% log_message('Finished estimating shifts');
% 
% % Reconstruct downsampled volume with no CTF correction
% n=size(projs,1);
% [ v1, ~, ~ ,~, ~, ~] = recon3d_firm( projs,...
%     rotations,-est_shifts, 1e-6, 100, zeros(n,n,n));
% ii1=norm(imag(v1(:)))/norm(v1(:));
% log_message('Relative norm of imaginary components = %e\n',ii1);
% v1=real(v1);
% WriteMRC(v1,1,outvol);
% 
% log_silent(currentsilentmode);