% A call to impl with only mandatory variables
n_symm = 5;
% mrc_stack_file = '/scratch/yoel/denoised/10089/denoised_group1.mrcs';
mrc_stack_file = '/home/gabip/matlabProjects/aspire/aspire/development/abinitio/C5/ml_cn/mrc_test_small.mrc';
recon_folder = '/home/gabip/matlabProjects/aspire/aspire/development/abinitio/C5/ml_cn/results_tmp';
recon_mrc_fname = fullfile(recon_folder,'out.mrc');
log_fname = fullfile(recon_folder,'log.txt');
open_log(log_fname);

% cryo_abinitio_Cn_ml_execute(n_symm,mrc_stack_file,recon_mrc_fname);

cache_file_name = '/home/gabip/matlabProjects/aspire/aspire/development/abinitio/C5/ml_cn/ml_cn_cache_points1000_ntheta360_res1.mat';
cryo_abinitio_Cn_ml_execute(n_symm,mrc_stack_file,recon_mrc_fname,cache_file_name);

% %% A call to impl with all possible variables
% cache_file_name = '/home/gabip/matlabProjects/aspire/aspire/development/abinitio/C5/ml_cn/ml_cn_cache_points1000_ntheta360_res1.mat';
% recon_mat_fname = '/home/gabip/matlabProjects/aspire/aspire/development/abinitio/C5/ml_cn/results_tmp/C5_10089_c5_nims100.mat';
% do_downsample = true;
% downsample_size = 65;
% n_r_perc = 50;
% max_shift_perc = 15;
% shift_step = 0.5;
% mask_radius_perc = 70;
% inplane_rot_res = 1;
% is_conjugate_with_vii = true;
% 
% cryo_abinitio_Cn_ml_execute(n_symm,mrc_stack_file,recon_mrc_fname,cache_file_name,recon_mat_fname,...
%     do_downsample,downsample_size,n_r_perc,max_shift_perc,shift_step,mask_radius_perc,inplane_rot_res,is_conjugate_with_vii);

close_log();