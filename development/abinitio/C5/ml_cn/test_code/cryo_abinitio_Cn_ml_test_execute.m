function [err_in_degrees,mse] = cryo_abinitio_Cn_ml_test_execute(n_symm,recon_mrc_fname,cache_file_name,snr,n_images,...
    n_r_perc,max_shift_perc,shift_step,mask_radius_perc,do_handle_equators,inplane_rot_res)

% open_log(fullfile(workflow.info.working_dir, workflow.info.logfile));
proj_size = 65;
max_shift = ceil(proj_size*max_shift_perc/100);
[projs,refq,~,~] = generate_cn_images(n_symm,n_images,snr,proj_size,'C1_Eytan',max_shift,shift_step);

[projs,refq] = remove_eq_images(projs,refq);
n_images = size(refq,2);

log_message('loading line indeces cache %s.\n Please be patient...',cache_file_name);
load(cache_file_name);
log_message('done loading indeces cache');

if snr <= 1
    mask_radius = proj_size*mask_radius_perc/100;
    log_message('Masking images using mask-radius=%d',mask_radius);
    masked_projs = mask_fuzzy(projs,mask_radius);
    
else
    masked_projs = projs;
    log_message('SNR=%.2e is greater than 1. Not performing mask', snr);
end

precision = 'double';
n_r = ceil(proj_size*n_r_perc/100);
[npf,~] = cryo_pft(masked_projs,n_r,n_theta,precision);

if snr <= 1
    log_message('Guass filtering the images');
    npf = gaussian_filter_imgs(npf);
else
    log_message('SNR=%.2e is greater than 1. Not performing gauss filtering', snr);
end

log_message('computing self common-line indeces for all candidate viewing directions');
ciis = compute_scls_inds(Ris_tilde,n_symm,n_theta);

is_viz_cls = false;
[vijs,viis,~] = compute_third_row_outer_prod_both_cn(npf,ciis,cijs_inds,Ris_tilde,R_theta_ijs,n_symm,max_shift,shift_step,...
    do_handle_equators,refq,is_viz_cls);

vijs = mat2flat(vijs,n_images);
[vijs,viis,~,~] = global_sync_J(vijs,viis);
% 
vis  = estimate_third_rows_ml(vijs,viis);
rots = estimate_inplane_rotations(npf,vis,n_symm,inplane_rot_res,max_shift,shift_step);

[rot_alligned,err_in_degrees,mse] = analyze_results_ml(rots,n_symm,n_theta,refq);
%
log_message('Reconstructing abinitio volume');
estimatedVol = reconstruct_ml_cn(projs,rot_alligned,n_symm,n_r,n_theta,max_shift,shift_step);
save_vols(estimatedVol,recon_mrc_fname,n_symm);

close_log();

end