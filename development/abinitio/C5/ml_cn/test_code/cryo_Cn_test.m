% % %% using only required input variables
% n_symm = 4;
% n_images = 100;
% recon_folder = '/home/gabip/matlabProjects/aspire/aspire/development/abinitio/C5/ml_cn/test_code/results';
% mrc_fname = sprintf('out_c%dnims%d.mrc',n_symm,n_images);
% recon_mrc_fname = fullfile(recon_folder,mrc_fname);
% 
% [err_in_degrees,mse] = cryo_abinitio_Cn_ml_test_execute(n_symm,recon_mrc_fname);
% % 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% using only required input variables + cache file
% clear;
% n_symm = 11;
% n_images = 100;
% recon_folder = '/home/gabip/matlabProjects/aspire/aspire/development/abinitio/C5/ml_cn/test_code/results';
% mrc_fname = sprintf('out_c%dnims%d.mrc',n_symm,n_images);
% recon_mrc_fname = fullfile(recon_folder,mrc_fname);
% 
% cache_file_name = '/home/gabip/matlabProjects/aspire/aspire/development/abinitio/C5/ml_cn/ml_cn_cache_points1000_ntheta360_res1.mat';
% [err_in_degrees,mse] = cryo_abinitio_Cn_ml_test_execute(n_symm,recon_mrc_fname,cache_file_name);
% % 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %% using all variables
% clear;
% n_symm = 4;
% snr = 100000000000;
% n_images = 100;
% n_r_perc = 50;
% max_shift_perc = 15;
% shift_step = 0.5;
% mask_radius_perc = 70;
% do_handle_equators = false;
% inplane_rot_res = 1;
% 
% recon_folder = '/home/gabip/matlabProjects/aspire/aspire/development/abinitio/C5/ml_cn/test_code/results';
% mrc_fname = sprintf('out_c%dnims%dshift%d.mrc',n_symm,n_images,max_shift_perc);
% recon_mrc_fname = fullfile(recon_folder,mrc_fname);
% cache_file_name = '/home/gabip/matlabProjects/aspire/aspire/development/abinitio/C5/ml_cn/ml_cn_cache_points1000_ntheta360_res1.mat';
% 
% 
% [err_in_degrees,mse] = cryo_abinitio_Cn_ml_test_execute(n_symm,recon_mrc_fname,cache_file_name,snr,n_images,...
%     n_r_perc,max_shift_perc,shift_step,mask_radius_perc,do_handle_equators,inplane_rot_res);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% a for loop that iterates over differernt symmetry classes
clear;
n_images = 100;
n_r_perc = 50;
max_shift_perc = 15;
recon_folder = '/home/gabip/matlabProjects/aspire/aspire/development/abinitio/C5/ml_cn/test_code/results';
cache_file_name = '/home/gabip/matlabProjects/aspire/aspire/development/abinitio/C5/ml_cn/ml_cn_cache_points1000_ntheta360_res1.mat';
shift_step = 0.5;
mask_radius_perc = 70;
% do_handle_equators = false;
inplane_rot_res = 1;

symms_s = 4;
snr_s = [100000000000];

res_med_err_in_degrees = zeros(numel(snr_s),numel(symms_s));
res_mse = zeros(numel(snr_s),numel(symms_s));
for i_symm = 1:numel(symms_s)
    for i_snr = 1:numel(snr_s)
        
        n_symm = symms_s(i_symm);
        snr = snr_s(i_snr);

        mrc_fname = sprintf('out_c%dnims%dshift%d.mrc',n_symm,n_images,max_shift_perc);
        recon_mrc_fname = fullfile(recon_folder,mrc_fname);
        
        [err_in_degrees,mse] = cryo_abinitio_Cn_ml_test_execute(n_symm,recon_mrc_fname,cache_file_name,snr,n_images,...
            n_r_perc,max_shift_perc,shift_step,mask_radius_perc,inplane_rot_res);
        
        res_med_err_in_degrees(i_snr,i_symm) = median(err_in_degrees);
        res_mse(i_snr,i_symm) = mse;
        
    end
end