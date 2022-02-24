function cryo_abinitio_TO(sym,instack,outvol,cache_file_name,matfname,...
          n_theta,n_r_perc,max_shift_perc,shift_step,mask_radius_perc,refq)
% cryo_abinitio_TO  abinitio reconsturction of a T or O symmetric molecule
%
% Parameters
%     sym             'T' for tetrahedral symmetry. 'O' for octahedral symmetry.
%   instack           Name of MRC file containing the projections (or class
%                     averages) from which to estimate an abinitio model.
%   outvol            Name of MRC file into which to save the reconstructed volume.
% cache_file_name     (Optional) The mat file name containing all candidate rotation
%                     matrices, common lines indices, self common lines indices. 
%                     If not supplied cache will be created.  
%   matfname          (Optional) Name of MAT file in which intermediate outcomes of the
%                     reconstruction algorithm are save (such as estimated
%                     rotations and shifts). Used to debugging and detailed
%                     analysis of the results.
%   n_theta           (Optional) Angular resolution for common lines detection.
%                     Default 360.
%   n_r_perc          (Optional) Radial resolution for common line detection as a
%                     percentage of image size. Default is half the width of the images.
%   max_shift_perc    (Optional) Maximal 1d shift (in pixels) to search between
%                     common-lines. Default is 15% of image width of the images.
%   shift_step        (Optional) Resolution of shift estimation in pixels. Note
%                     that shift_step can be any positive real number. Default:0.5.
%   mask_radius_perc  (Optional) The size of mask applied to each image as a
%                     percentage of image's size
%   refq              (Optional) A 4xnImages array holding the ground truth quaternions of all rotation matrices 
%
%  Adi Shasha, June 2021. 

if sym~='T' && sym~='O'
    error('First argument must be ''T'' or ''O'' (tetrahedral or octahedral symmetry)');
end

[folder_recon_mrc_fname, ~, ~] = fileparts(outvol);
if ~isempty(folder_recon_mrc_fname)  && exist(folder_recon_mrc_fname,'file') ~= 7
    error('folder %s does not exist. Please create it first.\n', folder_recon_mrc_fname);
end

if ~exist('cache_file_name','var') || isempty(cache_file_name) || ~exist(cache_file_name, 'file')
    log_message('Cache file not supplied.');
    [folder, ~, ~] = fileparts(outvol);
    if isempty(folder)
        folder = pwd;
    end
    cache_dir_full_path = folder;
    log_message('Creating cache file under folder: %s',cache_dir_full_path);
    cache_file_name = cryo_TO_create_cache(sym);
end

if exist('refq','var') && ~isempty(refq)
    refq_is_given = true;
else
    refq_is_given = false;
    refq = [];
end

if exist('matfname','var')
    do_save_res_to_mat = true;
else
    do_save_res_to_mat = false;
    matfname = '';
end

if ~exist('n_theta','var')
    n_theta = 360;
end

if ~exist('n_r_perc','var')
    n_r_perc = 50;
end

if ~exist('max_shift_perc','var')
    max_shift_perc = 15;
end
if ~exist('shift_step','var')
    shift_step = 0.5;
end

if ~exist('mask_radius_perc','var')
    mask_radius_perc = 50;
    
end

initstate; 


% Load projections
log_message('Loading mrc image stack file: %s.', instack);
projs = ReadMRC(instack);
projs = projs(1:end-1,1:end-1,:);
n_images = size(projs,3);
log_message('Done loading mrc image stack file');
log_message('projections loaded. Using %d projections of size %d x %d',n_images,size(projs,1),size(projs,2));


% figure; viewstack(projs,5,5);
mask_radius = ceil(size(projs,1)*mask_radius_perc/100);
if mask_radius > 0
    log_message('Masking projections. Masking radius is %d pixels',mask_radius);
    masked_projs = mask_fuzzy(projs, mask_radius);
else
    masked_projs = projs;
end

if do_save_res_to_mat
    log_message('Saving masked_projs under: %s', matfname);
    save(matfname,'masked_projs');
end

n_r = ceil(size(projs,1)*n_r_perc/100);
log_message('Computing the polar Fourier transform of projections');
[npf,~] = cryo_pft(masked_projs,n_r,n_theta,'double');
log_message('gassian filtering images');
npf = gaussian_filter_imgs(npf);
pf_norm = cryo_raynormalize(npf);

max_shift = ceil(size(projs,1)*max_shift_perc/100);
log_message('Maximum shift is %d pixels',max_shift);
log_message('Shift step is %d pixels',shift_step);

%% Step 3: Computing the relative rotations
log_message('Computing all relative rotations');
est_rel_rots = estimate_relative_rotations_with_shifts(pf_norm, max_shift, shift_step, cache_file_name);

if do_save_res_to_mat
    log_message('Saving relative rotations under: %s', matfname);
    save(matfname,'est_rel_rots','-append');
end


%% Step 4: Handedness synchronization
log_message('Handedness synchronization');
u_G = handedness_synchronization_TO(sym, est_rel_rots, cache_file_name);

if do_save_res_to_mat
    log_message('Saving handedness synchronization vector: %s', matfname);
    save(matfname,'u_G','-append');
end

%% Step 5: Rotation estimation
log_message('Estimating rotations');
rots = estimate_rotations_synchronization(est_rel_rots, u_G, cache_file_name);
rots_t = permute(rots,[2 1 3]);

if do_save_res_to_mat
    log_message('Saving rotations under: %s', matfname);
    save(matfname,'rots','-append');
end

%% Step 7: Results Analysis
if refq_is_given
    % the output rots are alligned with the ground-truth rotations
    [rots,err_in_degrees,mse] = analyze_results_TO(rots,sym,n_theta,refq);
end

%% Step 8: Reconstructing volume
log_message('Reconstructing abinitio volume');
[estimatedVol,est_rotations,est_shifts] = cryo_reconstruct_TO(sym,projs,rots_t,n_r,n_theta,max_shift,shift_step);

%% Step 9: Saving volume
WriteMRC(estimatedVol,1,outvol);

if do_save_res_to_mat
    log_message('Saving estimated volume, rotations and shifts under: %s', matfname);
    save(matfname,'estimatedVol','est_rotations','est_shifts','-append');
end

% close_log();

end

