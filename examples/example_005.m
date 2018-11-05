% This script illustrates the basic ab initio reconstruction functionality of
% the ASPIRE toolbox on Cn symmetric simulated data. 

% %% using all variables
clear;
n_symm = 3; % order of symmetry in the simulation n > 2
n_theta = 360; 
snr = 100000000000;
n_images = 100; % number of images
n_r_perc = 50; % percent of image size
max_shift_perc = 10; % percent of images size
shift_step = 0.5;
mask_radius_perc = 50;
inplane_rot_res = 1; 
verbose = 1;

recon_folder = '.';
cache_file_name = 'cn_cache_points1000_ntheta360_res1.mat';
mrc_fname = sprintf('out_c%dnims%dshift%d.mrc',n_symm,n_images,max_shift_perc);
recon_mrc_fname = fullfile(recon_folder,mrc_fname);

[err_in_degrees,mse] = cryo_abinitio_cn_test_execute(cache_file_name,n_symm,n_theta,recon_mrc_fname,snr,n_images,n_r_perc,...
    max_shift_perc,shift_step,mask_radius_perc,inplane_rot_res,verbose);

