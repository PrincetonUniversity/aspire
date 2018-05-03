n_symm = 4;
cache_file_name = '/home/gabip/matlabProjects/aspire/aspire/development/abinitio/C5/ml_cn/ml_cn_cache_points1000_ntheta360_res1.mat';
snr = 100000000000;
n_images = 100;
n_r_perc = 50;
max_shift_perc = 15;
shift_step = 0.5;
mask_radius_perc = 70;
do_handle_equators = false;
inplane_rot_res = 1;
recon_folder = '/home/gabip/matlabProjects/aspire/aspire/development/abinitio/C5/ml_cn/test_code/results';
mrc_fname = sprintf('out_c%dnims%dshit%d',n_symm,n_images,max_shift_perc);
recon_mrc_fname = fullfile(recon_folder,mrc_fname);
[err_in_degrees,mse] = cryo_abinitio_Cn_ml_test_execute(n_symm,recon_mrc_fname,cache_file_name,snr,n_images,...
    n_r_perc,max_shift_perc,shift_step,mask_radius_perc,do_handle_equators,inplane_rot_res);