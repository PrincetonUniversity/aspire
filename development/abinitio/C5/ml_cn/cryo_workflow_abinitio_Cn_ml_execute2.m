function cryo_workflow_abinitio_Cn_ml_execute2(n_symm,mrc_stack_file,recon_mrc_fname,cache_file_name,recon_mat_fname,...
    do_downsample,downsample_size,n_r_perc,max_shift_perc,shift_step,mask_radius_perc,do_handle_equators,inplane_rot_res)

if ~exist('cache_file_name','var')
    n_Points_sphere = 100;
    n_theta = 360;
    inplane_rot_res = 1;
    [folder, ~, ~] = fileparts(recon_mrc_fname);
    cache_dir_full_path = folder;
    log_message('Creating cache file under folder: %s',cache_dir_full_path);
    log_message('#points on sphere=%d, n_theta=%d, inplane_rot_res=%d',n_Points_sphere,n_theta,inplane_rot_res);
    cache_file_name  = cryo_Cn_ml_create_cache_mat(cache_dir_full_path,n_Points_sphere,n_theta,inplane_rot_res);
end

if ~exist('recon_mat_fname','var')
    do_save_res_to_mat = false;
else
    do_save_res_to_mat = true;
end

if ~exist('downsample_size','var')
    if exist('do_downsample','var')
        error('must specify downsample size');
    end
end

if ~exist('do_downsample','var')
    do_downsample = false;
end

if ~exist('n_r_perc','var')
    n_r_perc = 50;
end
if ~exist('max_shift_perc','var')
    max_shift_perc = 0;
end
if ~exist('shift_step','var')
    shift_step = 0.5;
end
if ~exist('mask_radius_perc','var')
    mask_radius_perc = 70;
    
end
if ~exist('do_handle_equators','var')
    do_handle_equators = false;
end

if ~exist('inplane_rot_res','var')
    inplane_rot_res = 1;
end
initstate; 

log_message('symmetry class is C%d',n_symm);

% Load projections
log_message('Loading mrc image stack file:%s. Plese be patient...', mrc_stack_file);
% log_message('Loading %d images, starting from image index %d',nImages,first_image_ind);
projs = ReadMRC(mrc_stack_file,1,100);
nImages = size(projs,3);
log_message('done loading mrc image stack file');
log_message('projections loaded. Using %d projections of size %d x %d',nImages,size(projs,1),size(projs,2));

if do_downsample
    projs = cryo_downsample(projs,downsample_size,1);
    log_message('Downsampled projections to be of size %d x %d',size(projs,1),size(projs,2));
end

log_message('loading line indeces cache %s.\n Please be patient...',cache_file_name);
load(cache_file_name);
log_message('done loading indeces cache');

log_message('computing self common-line indeces for all candidate viewing directions');
ciis = compute_scls_inds(Ris_tilde,n_symm,n_theta);

% figure; viewstack(projs,5,5);
mask_radius = ceil(size(projs,1)*mask_radius_perc/100);
if mask_radius > 0
    log_message('Masking projections. Masking radius is %d pixels',mask_radius);
    masked_projs = mask_fuzzy(projs, mask_radius);
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
[vijs,viis] = compute_third_row_outer_prod_both_cn(npf,ciis,cijs_inds,Ris_tilde,R_theta_ijs,...
    n_symm,max_shift,shift_step,do_handle_equators);

if do_save_res_to_mat
    log_message('Saving third rows outer prods under: %s', recon_mat_fname);
    save(recon_mat_fname,'vijs','viis');
end
% [vijs,viis,max_corrs_stats] = compute_third_row_outer_prod_both(npf,ciis,cijs,Ris_tilde,R_theta_ijs,max_shift,shift_step,is_handle_equators);

vijs = mat2flat(vijs,nImages);
[vijs,viis] = global_sync_J(vijs,viis);
if do_save_res_to_mat
    log_message('Saving third rows outer prods under: %s', recon_mat_fname);
    save(recon_mat_fname,'vijs','viis');
end
% 
vis  = estimate_third_rows_ml(vijs,viis);
if do_save_res_to_mat
    log_message('Saving third rows under: %s', recon_mat_fname);
    save(recon_mat_fname,'vis','-append');
end

rots = estimate_inplane_rotations(npf,vis,n_symm,inplane_rot_res,max_shift,shift_step);
if do_save_res_to_mat
    log_message('Saving estimated rotation matrices under: %s', recon_mat_fname);
    save(recon_mat_fname,'rots','-append');
end

log_message('Reconstructing abinitio volume');
estimatedVol = reconstruct_ml_cn(masked_projs,rots,n_symm,n_r,n_theta,max_shift,shift_step);
save_vols(estimatedVol,recon_mrc_fname,n_symm);

% close_log();

end

