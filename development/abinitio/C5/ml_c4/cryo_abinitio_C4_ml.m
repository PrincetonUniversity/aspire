function cryo_abinitio_C4_ml(instack,outvol,outparams,...
    nImages,n_theta,n_r,max_shift,shift_step)

nPoints_sphere  = 1000;
is_handle_equators = false;
inplane_rot_res = 1;

initstate; 
open_log(0);

% Check input and set default parameters
if ~exist('n_theta','var')
    n_theta = 360;
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

projs = projs(:,:,floor(linspace(1,size(projs,3),nImages)));
nImages = size(projs,3);

log_message('projections loaded. Using %d projections of size %d x %d',nImages,size(projs,1),size(projs,2));
if size(projs,1)~=size(projs,2)
    error('Input images must be square');
end

figure; viewstack(projs,5,5);

% mask_radius = round(size(projs,1)*0.45);
mask_radius = 30;
log_message('Masking projections. Masking radius is %d pixels',mask_radius);
masked_projs = mask_fuzzy(projs,mask_radius);
    
if ~n_r_given
    n_r = ceil(size(masked_projs,1)*0.5);
end
[npf,~] = cryo_pft(masked_projs,n_r,n_theta,'double');

npf = gaussian_filter_imgs(npf);

file_cache_name  = sprintf('ml_c4_cached_inds_%d.mat',nPoints_sphere);
load(file_cache_name);

% if is_use_gt_in_cands || is_save_inds_to_cache
%     [Ris_tilde,R_theta_ijs] = generate_cand_rots(nPoints_sphere,inplane_rot_res,is_use_gt_in_cands,refq);
%     [ciis,cijs]             = compute_cls_inds(Ris_tilde,R_theta_ijs,n_theta,is_save_inds_to_cache,file_cache_name);
% else
%     load(file_cache_name);
% end
if ~max_shift_given
    max_shift = ceil(size(projs,1)*0.15); % max_shift is 15% of the image size
%     max_shift = 3;
end

[vijs,viis,max_corrs_stats] = compute_third_row_outer_prod(npf,ciis,cijs,Ris_tilde,...
    R_theta_ijs,max_shift,shift_step,is_handle_equators);

[vijs,viis,sign_ij_J,sign_ii_J] = global_sync_J(vijs,viis);
% 
vis  = estimate_third_rows_ml(vijs,viis);
rots = estimate_inplane_rotations4(npf,vis,inplane_rot_res,max_shift,shift_step);

% [rot_alligned,err_in_degrees,mse] = analyze_results_ml(rots,n_theta,refq);

rot_alligned = rots;

estimatedVol = reconstruct_ml_c4(projs,rot_alligned,n_r,n_theta,max_shift,shift_step);
WriteMRC(estimatedVol,1,outvol);