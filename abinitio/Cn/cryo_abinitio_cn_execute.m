function cryo_abinitio_cn_execute(cache_file_name,n_symm,mrc_stack_file,recon_mrc_fname,recon_mat_fname,...
    verbose,n_theta,n_r_perc,max_shift_perc,shift_step,mask_radius_perc,inplane_rot_res)

[folder_recon_mrc_fname, ~, ~] = fileparts(recon_mrc_fname);
if ~isempty(folder_recon_mrc_fname)  && exist(folder_recon_mrc_fname,'file') ~= 7
    error('folder %s does not exist. Please create it first.\n', folder_recon_mrc_fname);
end


if ~exist('cache_file_name','var') || isempty(cache_file_name) || ~exist(cache_file_name, 'file')
    log_message('Cache file not supplied.');
    n_Points_sphere = 1000;
    n_theta = 360;
    inplane_rot_res = 1;
    [folder, ~, ~] = fileparts(recon_mrc_fname);
    cache_dir_full_path = folder;
    log_message('Creating cache file under folder: %s',cache_dir_full_path);
    log_message('#points on sphere=%d, n_theta=%d, inplane_rot_res=%d',n_Points_sphere,n_theta,inplane_rot_res);
    cache_file_name  = cryo_cn_create_cache(cache_dir_full_path,n_Points_sphere,n_theta,inplane_rot_res);
end


if ~exist('verbose','var')
    verbose = 0;
    log_message('verbose was not provided. Setting verbose=0\n');
end


if exist('recon_mat_fname','var')
    do_save_res_to_mat = true;
else
    do_save_res_to_mat = false;
    recon_mat_fname = '';
end

% if ~exist('downsample_size','var')
%     if exist('do_downsample','var')
%         error('must specify downsample size');
%     end
% end
% 
% if ~exist('do_downsample','var')
%     do_downsample = false;
% end

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


if(n_symm==3 || n_symm==4) % an empirical observation. But better check the alternative if reconstruction is not good enough
    is_conjugate_with_vii = false;
else
    is_conjugate_with_vii = true;
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
log_message('Loading mrc image stack file:%s. Plese be patient...', mrc_stack_file);
% log_message('Loading %d images, starting from image index %d',nImages,first_image_ind);
projs = ReadMRC(mrc_stack_file);
% projs = projs(:,:,ceil(linspace(1,5000,1500)));
nImages = size(projs,3);
log_message('done loading mrc image stack file');
log_message('projections loaded. Using %d projections of size %d x %d',nImages,size(projs,1),size(projs,2));

% if do_downsample
%     projs = cryo_downsample(projs,downsample_size,1);
%     log_message('Downsampled projections to be of size %d x %d',size(projs,1),size(projs,2));
% end

% figure; viewstack(projs,5,5);
mask_radius = ceil(size(projs,1)*mask_radius_perc/100);
if mask_radius > 0
    log_message('Masking projections. Masking radius is %d pixels',mask_radius);
    masked_projs = mask_fuzzy(projs, mask_radius);
else
    masked_projs = projs;
end

if do_save_res_to_mat
    log_message('Saving masked_projs under: %s', recon_mat_fname);
    save(recon_mat_fname,'masked_projs');
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


if(n_symm==3 || n_symm==4)
%     max_shift_1d  = ceil(2*sqrt(2)*max_shift); % TODO: is this needed? if so, where?
    is_remove_non_rank1 = true;
    non_rank1_remov_percent = 0.25;
    [vijs,viis,npf,masked_projs] = compute_third_row_outer_prod_c34(n_symm,npf,max_shift,shift_step,recon_mat_fname,...
        masked_projs,verbose,is_remove_non_rank1,non_rank1_remov_percent);
else
    [vijs,viis] = compute_third_row_outer_prod_cn(npf,n_symm,max_shift,shift_step,cache_file_name,verbose);
end

if do_save_res_to_mat
    log_message('Saving third rows outer prods under: %s', recon_mat_fname);
    save(recon_mat_fname,'vijs','viis','-append');
end

[vijs,viis] = global_sync_J(vijs,viis);
if do_save_res_to_mat
    log_message('Saving third rows outer prods under: %s', recon_mat_fname);
    save(recon_mat_fname,'vijs','viis','-append');
end
% 
vis  = estimate_third_rows(vijs,viis,is_conjugate_with_vii);
if do_save_res_to_mat
    log_message('Saving third rows under: %s', recon_mat_fname);
    save(recon_mat_fname,'vis','-append');
end

rots = estimate_inplane_rotations(npf,vis,n_symm,inplane_rot_res,max_shift,shift_step,verbose);
if do_save_res_to_mat
    log_message('Saving estimated rotation matrices under: %s', recon_mat_fname);
    save(recon_mat_fname,'rots','-append');
end

log_message('Reconstructing abinitio volume');
estimatedVol = reconstruct_cn(masked_projs,rots,n_symm,n_r,n_theta,max_shift,shift_step);
save_vols(estimatedVol,recon_mrc_fname,n_symm);

if do_save_res_to_mat
    log_message('Saving estimated rotations and shoifts under: %s', recon_mat_fname);
    save(recon_mat_fname,'rotations','est_shifts','-append');
end


% close_log();

end
