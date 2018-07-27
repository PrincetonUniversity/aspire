function [vijs,viis,npf,projs,refq] = estimate_third_rows_outer_prod_c3_c4(n_symm,npf,max_shift,shift_step,recon_mat_fname,...
    projs,is_remove_non_rank1,non_rank1_remov_percent,refq)

if exist('refq','var') && ~isempty(refq)
    is_simulation = true;
else
    refq = [];
    is_simulation = false;
end

if exist('recon_mat_fname','var') && ~isempty(recon_mat_fname)
    do_save_res_to_mat = true;
else
    do_save_res_to_mat = false;
end

[~,n_theta,n_images] = size(npf);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 1  : detect a single pair of common-lines between each pair of images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
log_message('Detecting common-lines');
clmatrix = cryo_clmatrix(npf,n_images,1,max_shift,shift_step); 
if do_save_res_to_mat
    log_message('Saving clmatrix to %s',recon_mat_fname);
    save(recon_mat_fname,'clmatrix','-append');
end

if is_simulation
    cl_detection_rate_c3_c4(n_symm,clmatrix,n_theta,refq);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 2  : detect self-common-lines in each image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if n_symm == 3
    is_handle_equator_ims = false;
else
    is_handle_equator_ims = true;
end
sclmatrix = cryo_self_clmatrix_gpu_c3_c4(n_symm,npf,max_shift,shift_step,is_handle_equator_ims);
if do_save_res_to_mat
    log_message('Saving sclmatrix to %s',recon_mat_fname);
    save(recon_mat_fname,'sclmatrix','-append');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 3  : calculate self-relative-rotations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Riis = estimate_all_Riis_c3_c4(n_symm,sclmatrix,n_theta,refq);
if do_save_res_to_mat
    log_message('Saving Riis to %s',recon_mat_fname);
    save(recon_mat_fname,'Riis','-append');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 4  : calculate relative-rotations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rijs = cryo_estimate_all_Rijs_c3_c4(n_symm,clmatrix,n_theta,refq);
if do_save_res_to_mat
    log_message('Saving Rijs to %s',recon_mat_fname);
    save(recon_mat_fname,'Rijs','-append');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 5  : inner J-synchronization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[vijs,viis,im_inds_to_remove,pairwise_inds_to_remove,...
    npf,projs,refq] = local_sync_J_c3_c4(n_symm,Rijs,Riis,npf,projs,is_remove_non_rank1,non_rank1_remov_percent,refq);
if do_save_res_to_mat
    log_message('Saving npf to %s',recon_mat_fname);
    log_message('Saving vijs to %s',recon_mat_fname);
    log_message('Saving viis to %s',recon_mat_fname);
    log_message('Saving im_inds_to_remove to %s',recon_mat_fname);
    log_message('Saving pairwise_inds_to_remove to %s',recon_mat_fname);
    save(recon_mat_fname,'npf','vijs','viis','im_inds_to_remove','pairwise_inds_to_remove','-append');
end

end