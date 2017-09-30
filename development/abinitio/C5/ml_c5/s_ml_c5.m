tic;
close all;
nImages = 100;
nPoints_sphere  = 1000;
is_use_gt_in_cands = false;
is_save_inds_to_cache = false;
n_r     = 65;  
n_theta = 360;
snr = 100000000000;
inplane_rot_res = 1;
c5_type = 'C1';
is_viz_cls = false;
max_shift  = 0;
shift_step = 0.5;

initstate; 
open_log(0)

[projs,refq] = generate_c5_images(nImages,snr,65,c5_type,max_shift,shift_step);
nImages = size(refq,2);

figure; viewstack(projs,5,5);

masked_projs = mask_fuzzy(projs,23);

[npf,~]      = cryo_pft(masked_projs,n_r,n_theta,'double');

npf = gaussian_filter_imgs(npf);

file_cache_name  = sprintf('ml_c5_cached_inds_%d.mat',nPoints_sphere);

if is_use_gt_in_cands || is_save_inds_to_cache
    [Ris_tilde,R_theta_ijs] = generate_cand_rots(nPoints_sphere,inplane_rot_res,is_use_gt_in_cands,refq);
    [ciis,cijs]             = compute_cls_inds(Ris_tilde,R_theta_ijs,n_theta,is_save_inds_to_cache,file_cache_name);
else
    load(file_cache_name);
end

[vijs,viis,max_corrs_stats] = compute_third_row_outer_prod(npf,ciis,cijs,Ris_tilde,R_theta_ijs,max_shift,shift_step,refq,is_viz_cls);

[vijs,viis,sign_ij_J,sign_ii_J] = global_sync_J(vijs,viis);
% 
vis  = estimate_third_rows_ml_c5(vijs,viis);
rots = estimate_inplane_rotations4(npf,vis,inplane_rot_res,max_shift,shift_step);

[rot_alligned,err_in_degrees,mse] = analyze_results_ml(rots,n_theta,refq);

toc

estimatedVol = reconstruct_ml_c5(projs,rot_alligned,n_r,n_theta,max_shift,shift_step);
WriteMRC(estimatedVol,1,'ml_cn.mrc');