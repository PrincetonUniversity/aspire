close all;
n_symm = 11; % symmetry class
nImages = 100;
nPoints_sphere  = 1000;
is_use_gt_in_cands = false;
is_save_inds_to_cache = false;
n_r     = 65;  
n_theta = 360;
snr = 1/8;
is_handle_equators = false;
inplane_rot_res_cache = 1;
inplane_rot_res = 1;
is_viz_cls = false;
max_shift  = 0;
shift_step = 0.5;

initstate;
open_log(0);

if is_handle_equators == true
    assert(mod(n_symm,2) == 0);
end

[projs,refq,~,clean_images] = generate_cn_images(n_symm,nImages,snr,65,'C1_Eytan',max_shift,shift_step);

[projs,refq] = remove_eq_images(projs,refq);

nImages = size(refq,2);

figure; viewstack(projs,5,5);

if snr <=1
    mask_radius = round(size(projs,1)*0.45);
    masked_projs = mask_fuzzy(projs,mask_radius);
else
    masked_projs = projs;
end

[npf,~]      = cryo_pft(masked_projs,n_r,n_theta,'double');

if snr <=1
    npf = gaussian_filter_imgs(npf);
end

file_cache_name  = sprintf('ml_c%d_cached_inds_%d.mat',n_symm,nPoints_sphere);
% file_cache_name  = sprintf('ml_c12_cached_inds_no_eq_%d.mat',nPoints_sphere);
% file_cache_name  = sprintf('ml_c12_cached_inds_%d_res%d.mat',nPoints_sphere,inplane_rot_res_cache);

if is_use_gt_in_cands || is_save_inds_to_cache
    [Ris_tilde,R_theta_ijs] = generate_cand_rots(nPoints_sphere,n_symm,inplane_rot_res_cache,is_use_gt_in_cands,refq);
    [ciis,cijs]             = compute_cls_inds(Ris_tilde,R_theta_ijs,n_symm,n_theta,is_save_inds_to_cache,file_cache_name);
else
    log_message('loading line indeces cache');
    load(file_cache_name);
end

[vijs,viis,max_corrs_stats] = compute_third_row_outer_prod_both2(npf,ciis,cijs,Ris_tilde,R_theta_ijs,n_symm,max_shift,shift_step,...
    is_handle_equators,refq,is_viz_cls);

vijs = mat2flat(vijs,nImages);
[vijs,viis,sign_ij_J,sign_ii_J] = global_sync_J(vijs,viis);
% 
vis  = estimate_third_rows_ml(vijs,viis);
rots = estimate_inplane_rotations(npf,vis,n_symm,inplane_rot_res,max_shift,shift_step);

[rot_alligned,err_in_degrees,mse] = analyze_results_ml(rots,n_symm,n_theta,refq);

% 
estimatedVol = reconstruct_ml_cn(projs,rot_alligned,n_symm,n_r,n_theta,max_shift,shift_step);
err = assertVolumeIsCn(estimatedVol,n_symm);
WriteMRC(estimatedVol,1,'ml_cn_from_noisy.mrc');

% 
% estimatedVol_clean = reconstruct_ml_cn(clean_images,rot_alligned,n_symm,n_r,n_theta,max_shift,shift_step);
% WriteMRC(estimatedVol_clean,1,'ml_cn_from_clean.mrc');