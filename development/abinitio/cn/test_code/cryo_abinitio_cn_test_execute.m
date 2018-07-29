function [err_in_degrees,mse] = cryo_abinitio_cn_test_execute(cache_file_name,n_symm,n_theta,recon_mrc_fname,snr,n_images,n_r_perc,...
    max_shift_perc,shift_step,mask_radius_perc,inplane_rot_res)

[folder_recon_mrc_fname, ~, ~] = fileparts(recon_mrc_fname);
if ~isempty(folder_recon_mrc_fname)  && exist(folder_recon_mrc_fname,'file') ~= 7
    error('Folder %s does not exist. Please create it first.\n', folder_recon_mrc_fname);
end

if ~exist('snr','var')
    snr = 10000000;
    log_message('snr not specified. Using %.2e\n',snr);
end

if ~exist('n_images','var')
    n_images = 100;
    log_message('n_images not specified. Using %d\n',n_images);
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
    mask_radius_perc = 70;
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
proj_size = 65;
max_shift = ceil(proj_size*max_shift_perc/100);

[projs,refq,~,~,vol_orig] = generate_cn_images(n_symm,n_images,snr,proj_size,'C1_Eytan',max_shift,shift_step);
% saving original volume to disk
[folder, ~, ~] = fileparts(recon_mrc_fname);
vol_orig_file_name = fullfile(folder,sprintf('vol_orig_c%d.mrc',n_symm));
log_message('saving original volume under %s',vol_orig_file_name);
WriteMRC(vol_orig,1,vol_orig_file_name);

if n_symm > 4
    % currently we don't know how to handle equators for n_symm > 4
    [projs,refq] = remove_eq_images(projs,refq);
end

if snr <= 1
    mask_radius = proj_size*mask_radius_perc/100;
    log_message('Masking images using mask-radius=%d',mask_radius);
    masked_projs = mask_fuzzy(projs,mask_radius);
    
else
    masked_projs = projs;
    log_message('SNR=%.2e is greater than 1. Not performing mask', snr);
end

n_r = ceil(proj_size*n_r_perc/100);
[npf,~] = cryo_pft(masked_projs,n_r,n_theta,'single');

if snr <= 1
    log_message('Guass filtering the images');
    npf = gaussian_filter_imgs(npf);
else
    log_message('SNR=%.2e is greater than 1. Not performing gauss filtering', snr);
end

if(n_symm==3 || n_symm==4)
    max_shift_1d  = ceil(2*sqrt(2)*max_shift); % TODO: is this needed? if so, where?
    is_remove_non_rank1 = true;
    non_rank1_remov_percent = 0.25;
    [vijs,viis,npf,projs,refq] = compute_third_row_outer_prod_c34(n_symm,npf,max_shift_1d,shift_step,'',...
        projs,is_remove_non_rank1,non_rank1_remov_percent,refq);
else
    [vijs,viis,~] = compute_third_row_outer_prod_cn(npf,n_symm,max_shift,shift_step,cache_file_name,refq);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 7  : outer J-synchronization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[vijs,viis] = global_sync_J(vijs,viis);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 8  : third rows estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vis = estimate_third_rows(vijs,viis,is_conjugate_with_vii);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 9  : in-plane rotations angles estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rots =      estimate_inplane_rotations(npf,vis,n_symm,inplane_rot_res,max_shift,shift_step);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 10  : Results Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[rot_alligned,err_in_degrees,mse] = analyze_results_cn(rots,n_symm,n_theta,refq);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 11  : Reconstructing volume
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
estimatedVol = reconstruct_cn(projs,rot_alligned,n_symm,n_r,n_theta); % supply max_shift_2d or max_shift_1d ?

save_vols(estimatedVol,recon_mrc_fname,n_symm);

end