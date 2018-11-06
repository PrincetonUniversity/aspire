function cryo_abinitio_Cn(n_symm,instack,outvol,cache_file_name,matfname,...
    verbose,n_theta,n_r_perc,max_shift_perc,shift_step,mask_radius_perc,inplane_rot_res,refq)
% cryo_abinitio_Cn  abinitio reconsturction of a Cn symmetric molecule
%
% Parameters
%   n_symm      The symmetry order of the underlying molecule (n_symm > 2)
%   instack     Name of MRC file containing the projections (or class
%               averages) from which to estimate an abinitio model.
%   outvol      Name of MRC file into which to save the reconstructed volume.
%cache_file_name (Optional) The mat file name containing all candidate rotation
%               matrices. Only used for n_symm > 4. If not supplied cache will be created.  
%   matfname    (Optional) Name of MAT file in which intermediate outcomes of the
%               reconstruction algorithm are save (such as estimated
%               rotations and shifts). Used to debugging and detailed
%               analysis of the results.
%   verbose     (Optional) Level of detailed debug output. 0: none, 1:
%                detailed output (default:0)
%   ntheta      (Optional) Angular resolution for common lines detection.
%               Default 360.
%   n_r_perc   (Optional) Radial resolution for common line detection as a
%               percentage of image size. Default is half the width of the images.
%   max_shift   (Optional) Maximal 1d shift (in pixels) to search between
%               common-lines. Default is 15% of image width of the images.
%   shift_step  (Optional) Resolution of shift estimation in pixels. Note
%               that shift_step can be any positive real number. Default:0.5.
%   mask_radius_perc (Optional) The size of mask applied to each image as a
%                   percentage of image's size
%   inplane_rot_res (Optional) the resolution in angles to search for the
%                   inplane rotation angle of each rotation matrix (Default:1)
%   refq          (Optional) A 4xnImages array holding the ground truth quaternions of all rotation matrices 
%

[folder_recon_mrc_fname, ~, ~] = fileparts(outvol);
if ~isempty(folder_recon_mrc_fname)  && exist(folder_recon_mrc_fname,'file') ~= 7
    error('folder %s does not exist. Please create it first.\n', folder_recon_mrc_fname);
end

if n_symm > 4 % check that cache file exists, and create it otherwise
    if ~exist('cache_file_name','var') || isempty(cache_file_name) || ~exist(cache_file_name, 'file')
        log_message('Cache file not supplied.');
        n_Points_sphere = 1000;
        n_theta = 360;
        inplane_rot_res = 1;
        [folder, ~, ~] = fileparts(outvol);
        cache_dir_full_path = folder;
        log_message('Creating cache file under folder: %s',cache_dir_full_path);
        log_message('#points on sphere=%d, n_theta=%d, inplane_rot_res=%d',n_Points_sphere,n_theta,inplane_rot_res);
        cache_file_name  = cryo_cn_create_cache(cache_dir_full_path,n_Points_sphere,n_theta,inplane_rot_res);
    end
end

if exist('refq','var') && ~isempty(refq)
    refq_is_given = true;
else
    refq_is_given = false;
    refq = [];
end

if ~exist('verbose','var')
    verbose = 0;
    log_message('verbose was not provided. Setting verbose=0\n');
end


if exist('recon_mat_fname','var')
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

if ~exist('inplane_rot_res','var')
    inplane_rot_res = 1;
end

initstate; 

log_message('symmetry class is C%d',n_symm);


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
% Load projections
log_message('Loading mrc image stack file:%s. Plese be patient...', instack);
% log_message('Loading %d images, starting from image index %d',nImages,first_image_ind);
projs = ReadMRC(instack);
% projs = projs(:,:,ceil(linspace(1,5000,1500)));
nImages = size(projs,3);
log_message('done loading mrc image stack file');
log_message('projections loaded. Using %d projections of size %d x %d',nImages,size(projs,1),size(projs,2));


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

log_message('determining third rows outer product using maximum likelihood');
max_shift = ceil(size(projs,1)*max_shift_perc/100);
log_message('Maximum shift is %d pixels',max_shift);
log_message('Shift step is %d pixels',shift_step);

%% Step 3: Computing the relative viewing directions
if(n_symm==3 || n_symm==4)
    max_shift_1d  = ceil(2*sqrt(2)*max_shift); % TODO: is this needed? if so, where?
    is_remove_non_rank1 = true; % whether or not to remove images with the least number 
    %of induced rank-1 matrices of relative viewing directions 
    non_rank1_remov_percent = 0.25; % the percent of images to remove
    log_message('Computing all relative viewing directions for n=3,4');
    [vijs,viis,npf,masked_projs,refq] = compute_third_row_outer_prod_c34(n_symm,npf,max_shift_1d,shift_step,matfname,...
        masked_projs,verbose,is_remove_non_rank1,non_rank1_remov_percent,refq);
else
    log_message('Computing all relative viewing directions for n>4');
    [vijs,viis] = compute_third_row_outer_prod_cn(npf,n_symm,max_shift,shift_step,cache_file_name,verbose,refq);
end

if do_save_res_to_mat
    log_message('Saving third rows outer prods under: %s', matfname);
    save(matfname,'vijs','viis','-append');
end

%% Step 4: Handedness synchronization
log_message('Handedness synchronization');
[vijs,viis] = global_sync_J(vijs,viis);
if do_save_res_to_mat
    log_message('Saving third rows outer prods under: %s', matfname);
    save(matfname,'vijs','viis','-append');
end

%% Step 5: Viewing direction estimation (i.e., estimating third row of each rotation matrix)
log_message('Estimating the viewing direction (i.e. third row) of each image');
vis  = estimate_third_rows(vijs,viis,n_symm);
if do_save_res_to_mat
    log_message('Saving third rows under: %s', matfname);
    save(matfname,'vis','-append');
end

%% step 6: In-plane rotations angles estimation
log_message('In plane angles estimation');
rots = estimate_inplane_rotations(npf,vis,n_symm,inplane_rot_res,max_shift,shift_step,verbose);
if do_save_res_to_mat
    log_message('Saving estimated rotation matrices under: %s', matfname);
    save(matfname,'rots','-append');
end

%% Step 7: Results Analysis
if refq_is_given
    % the output rots are alligned with the ground-truth rotations
    [rots,err_in_degrees,mse] = analyze_results_cn(rots,n_symm,n_theta,refq);
end

%% Step 8: Reconstructing volume
log_message('Reconstructing abinitio volume');
estimatedVol = reconstruct_cn(masked_projs,rots,n_symm,n_r,n_theta,max_shift,shift_step); % supply max_shift_2d or max_shift_1d ?

%% Step 9: Saving volumes
save_vols(estimatedVol,outvol,n_symm);

if do_save_res_to_mat
    log_message('Saving estimated rotations and shoifts under: %s', matfname);
    save(matfname,'rotations','est_shifts','-append');
end

% close_log();

end