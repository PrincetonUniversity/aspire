function cryo_abinitio_C4(instack,outvol,outmat,max_shift_perc,shift_step,n_r_perc,mask_radius_perc,n_theta)

% CRYO_ABINITO_C4  abinitio reconsturction of a c4 symmetric molecule
%
% Parameters
%   instack     Name of MRC file containing the projections (or class
%               averages) from which to estimate an abinitio model.
%   outvol      Name of MRC file into which to save the reconstructed
%               volume.
%   outmat      (Optional) Name of MAT file in which intermediate outcomes of the
%               reconstruction algorithm are save (such as estimated
%               rotations and shifts). Used to debugging and detailed
%               analysis of the results.
%   ntheta      (Optional) Angular resolution for common lines detection.
%               Default 360. 
%   n_r_perc    (Optional) Radial resolution for common line detection as a
%               percentage of image size.
%               Default is half the width of the images.
%   max_shift_perc (Optional) Maximal 1d shift (as percentage of image size) 
%               to search between common-lines. Default is 15% of image width of the images.
%   shift_step  (Optional) Resolution of shift estimation in pixels. Note
%               that shift_step can be any positive real number. Default:0.5. 
%

% Check input and set default parameters
if ~exist('n_theta','var')
    n_theta = 360;
end

if ~exist('mask_radius_perc','var')
    mask_radius_perc = 70;
end


if ~exist('n_r_perc','var')
    n_r_perc = 50;
end

if ~exist('shift_step','var')
    shift_step = 0.5;
end

if ~exist('max_shift_perc','var')
    max_shift_perc = 15;
end

if ~exist('outmat','var')
    do_save_res_to_mat = false;
else
    do_save_res_to_mat = true;
end

if do_save_res_to_mat
    [folder_mat_fname, ~, ~] = fileparts(outmat);
    if ~isempty(folder_mat_fname)  && exist(folder_mat_fname,'file') ~= 7
        error('Folder %s does not exist. Please create it first.\n', folder_mat_fname);
    end
end

[folder_outvol_fname, ~, ~] = fileparts(outvol);
if ~isempty(folder_outvol_fname)  && exist(folder_outvol_fname,'file') ~= 7
    error('Folder %s does not exist. Please create it first.\n', folder_outvol_fname);
end


if do_save_res_to_mat
    log_message('Saving input variables to %s',outmat)
    save(outmat,'max_shift_perc','shift_step','n_r_perc','mask_radius_perc','n_theta');
end

% %% Load projections
% projs = ReadMRC(instack,1,5000);
% 
% if n_projs_given
%     if n_projs == -1
%         n_images = size(projs,3);
%     elseif n_projs <= 0 
%         error('n_projs must be either positive number, or -1 for using all images');
%     else
%         assert(n_projs <= size(projs,3));
%         n_images = n_projs;
%     end
% else % not provided as a parameter so use everything
%     n_images = size(projs,3);
% end
% 
% im_indeces = randperm(size(projs,3),n_images);
% save(outparams,'im_indeces');
% 
% projs = projs(:,:,im_indeces);
% assert(size(projs,3) == n_images);
% 
projs = ReadMRC(instack);
n_images = size(projs,3);
log_message('Read %d images',n_images);

% sz = size(ReadMRC(instack,1,1),1);
% projs = zeros(sz,sz,n_images);
% im_indeces  = randperm(40000,n_images);
% im_indeces  = sort(im_indeces);
% save(outparams,'im_indeces');
% projs = get_ims_from_stack(instack,im_indeces);

log_message('projections loaded. Using %d projections of size %d x %d',n_images,size(projs,1),size(projs,2));
if size(projs,1)~=size(projs,2)
    error('Input images must be square');
end

%% Mask projections
mask_radius = round(size(projs,1)*mask_radius_perc/100);
log_message('Masking projections. Masking radius is %d pixels',mask_radius);
[masked_projs,~] = mask_fuzzy(projs,mask_radius);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 1  : Computing polar Fourier transform of projections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_r = ceil(size(masked_projs,1)*n_r_perc/100);

[npf,~] = cryo_pft(masked_projs,n_r,n_theta,'single');

log_message('Applying Gaussian filter to images');
npf = gaussian_filter_imgs(npf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 2  : detect a single pair of common-lines between each pair of images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
max_shift = ceil(size(projs,1)*max_shift_perc/100);

log_message('Detecting common-lines');
clmatrix = cryo_clmatrix(npf,n_images,1,max_shift,shift_step); 
if do_save_res_to_mat
    log_message('Saving clmatrix to %s',outmat);
    save(outmat,'clmatrix','-append');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 3  : detect self-common-lines in each image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
is_handle_equator_ims = true;
sclmatrix = cryo_self_clmatrix_gpu(npf,max_shift,shift_step,is_handle_equator_ims);
if do_save_res_to_mat
    log_message('Saving sclmatrix to %s',outmat);
    save(outmat,'sclmatrix','-append');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 4  : calculate self-relative-rotations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Riis = estimate_all_Riis(sclmatrix,n_theta);
if do_save_res_to_mat
    log_message('Saving Riis to %s',outmat);
    save(outmat,'Riis','-append');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 5  : calculate relative-rotations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rijs = cryo_c4_estimate_all_Rijs(clmatrix,n_theta);
if do_save_res_to_mat
    log_message('Saving Rijs to %s',outmat);
    save(outmat,'Rijs','-append');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 6  : inner J-synchronization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[vijs,viis,im_inds_to_remove,pairwise_inds_to_remove,...
    npf,projs] = local_sync_J(Rijs,Riis,npf,projs);
if do_save_res_to_mat
    log_message('Saving npf to %s',outmat);
    log_message('Saving vijs to %s',outmat);
    log_message('Saving viis to %s',outmat);
    log_message('Saving im_inds_to_remove to %s',outmat);
    log_message('Saving pairwise_inds_to_remove to %s',outmat);
    save(outmat,'npf','vijs','viis','im_inds_to_remove','pairwise_inds_to_remove','-append');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 7  : outer J-synchronization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[vijs,viis] = global_sync_J(vijs,viis);
if do_save_res_to_mat
    log_message('Saving vijs to %s',outmat);
    log_message('Saving viis to %s',outmat);
    save(outmat,'vijs','viis','-append');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 8  : third rows estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vis  = estimate_third_rows(vijs,viis);
if do_save_res_to_mat
    log_message('Saving vis to %s',outmat);
    save(outmat,'vis','-append');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 9  : in-plane rotations angles estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rots = estimate_inplane_rotations4(npf,vis,1,max_shift,shift_step);
if do_save_res_to_mat
    log_message('Saving rots to %s',outmat);
    save(outmat,'rots','-append');
end

% estimatedVol = reconstruct_vol(projs,npf,rot_alligned,max_shift,shift_step);
estimatedVol = reconstruct(projs,rots,n_r,n_theta,max_shift,shift_step);   
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