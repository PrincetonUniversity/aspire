function cryo_workflow_abinitio_Cn_ml_execute(workflow_fname)

%% Validate workflow file
cryo_workflow_abinitio_Cn_ml_execute_validate(workflow_fname);

%% Read workflow file
tree = xmltree(workflow_fname);
workflow = convert(tree);

initstate; 
%% Execute preprocessing
open_log(fullfile(workflow.info.working_dir,workflow.info.logfile));

log_message('cryo_workflow_abinitio_Cn_ml_execute');
log_message('Loaded XML file %s',workflow_fname);

n_symm = str2double(workflow.algo.n_symm);
log_message('symmetry class is C%d',n_symm);
is_handle_equators = str2double(workflow.algo.do_handle_equators);

mrc_stack_file = workflow.info.mrc_stack_file;

first_image_ind = str2double(workflow.algo.first_image_ind);
last_image_ind = str2double(workflow.algo.last_image_ind);
nImages = last_image_ind - first_image_ind + 1;
% Load projections
log_message('Loading mrc image stack file:%s. Plese be patient...', mrc_stack_file);
log_message('Loading %d images, starting from image index %d',nImages,first_image_ind);
projs = ReadMRC(mrc_stack_file,first_image_ind,nImages);
log_message('done loading mrc image stack file');
log_message('projections loaded. Using %d projections of size %d x %d',nImages,size(projs,1),size(projs,2));

if str2double(workflow.algo.do_downsample)
    downsample_size = str2double(workflow.algo.downsample_size);
    projs = cryo_downsample(projs,downsample_size,1);
    log_message('Downsampled projections to be of size %d x %d',size(projs,1),size(projs,2));
end

log_message('loading line indeces cache %s.\n Please be patient...',workflow.cache.name);
load(workflow.cache.name);
log_message('done loading indeces cache');

log_message('computing self common-line indeces for all candidate viewing directions');
ciis = compute_scls_inds(Ris_tilde,n_symm,n_theta);

% figure; viewstack(projs,5,5);
mask_radius_perc = str2double(workflow.algo.mask_radius_perc);
mask_radius = ceil(size(projs,1)*mask_radius_perc/100);
if mask_radius > 0
    log_message('Masking projections. Masking radius is %d pixels',mask_radius);
    masked_projs = mask_fuzzy(projs, mask_radius);
end

n_r_perc = str2double(workflow.algo.n_r_perc);
n_r = ceil(size(projs,1)*n_r_perc/100);
log_message('Computing the polar Fourier transform of projections');
[npf,~] = cryo_pft(masked_projs,n_r,n_theta,'double');
log_message('gassian filtering images');
npf = gaussian_filter_imgs(npf);

log_message('determining third rows outer product using maximum likelihood');
max_shift_perc = str2double(workflow.algo.max_shift_perc);
max_shift = ceil(size(projs,1)*max_shift_perc/100);
log_message('Maximum shift is %d pixels',max_shift);
shift_step = str2double(workflow.algo.shift_step);
log_message('Shift step is %d pixels',shift_step);
[vijs,viis] = compute_third_row_outer_prod_both_cn(npf,ciis,cijs_inds,Ris_tilde,R_theta_ijs,...
    n_symm,max_shift,shift_step,is_handle_equators);

log_message('Saving third rows outer prods under: %s', workflow.algo.recon_mat_fname);
save(workflow.algo.recon_mat_fname,'vijs','viis');
% [vijs,viis,max_corrs_stats] = compute_third_row_outer_prod_both(npf,ciis,cijs,Ris_tilde,R_theta_ijs,max_shift,shift_step,is_handle_equators);

vijs = mat2flat(vijs,nImages);
[vijs,viis] = global_sync_J(vijs,viis);
log_message('Saving third rows outer prods under: %s', workflow.algo.recon_mat_fname);
save(workflow.algo.recon_mat_fname,'vijs','viis');
% 
vis  = estimate_third_rows_ml(vijs,viis);
log_message('Saving third rows under: %s', workflow.algo.recon_mat_fname);
save(workflow.algo.recon_mat_fname,'vis','-append');

inplane_rot_res = str2double(workflow.algo.inplane_rot_res);
rots = estimate_inplane_rotations(npf,vis,n_symm,inplane_rot_res,max_shift,shift_step);
log_message('Saving estimated rotation matrices under: %s', workflow.algo.recon_mat_fname);
save(workflow.algo.recon_mat_fname,'rots','-append');
log_message('Reconstructing abinitio volume');
estimatedVol = reconstruct_ml_cn(masked_projs,rots,n_symm,n_r,n_theta,max_shift,shift_step);

save_vols(estimatedVol,workflow.algo.recon_mrc_fname,n_symm);

close_log();

end

function cryo_workflow_abinitio_Cn_ml_execute_validate(workflow_fname)
% Read workflow file
tree = xmltree(workflow_fname);
workflow = convert(tree);

% Validate struct

assertfield(workflow,'cache','name');

is_valid = validate_inds_cache_file(workflow.cache.name);
if ~is_valid
    fprintf('Viewing directions inds cache %s is not valid\n',workflow.cache.name);
    fprintf('Please create new cache file using cryo_workflow_Cn_start');
    return;
end

assertfield(workflow,'info','name');
assertfield(workflow,'info','description');
assertfield(workflow,'info','mrc_stack_file');
assertfield(workflow,'info','created');
assertfield(workflow,'info','working_dir');
assertfield(workflow,'info','logfile');

assertfield(workflow,'algo','recon_mrc_fname');
assertfield(workflow,'algo','recon_mat_fname');
assertfield(workflow,'algo','n_symm');

assertfield(workflow,'algo','first_image_ind');
assertfield(workflow,'algo','last_image_ind');
assertfield(workflow,'algo','do_downsample');
assertfield(workflow,'algo','downsample_size');
% assertfield(workflow,'algo','n_images');
assertfield(workflow,'algo','n_r_perc');
assertfield(workflow,'algo','max_shift_perc');
assertfield(workflow,'algo','mask_radius_perc');
assertfield(workflow,'algo','shift_step');
assertfield(workflow,'algo','do_handle_equators');
assertfield(workflow,'algo','inplane_rot_res');

end