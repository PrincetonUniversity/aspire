% This script illustrates the basic ab initio reconstruction functionality of
% the ASPIRE toolbox on Cn symmetric simulated data. 

%% Parameters %%

n_symm = 3; % order of symmetry in the simulation (n > 2)
proj_size = 65; % the size in pixels of simulated projection images
n_theta = 360; % the angular resolution in nufft
snr = 100000000000; % the signal to noise ratio of the simulated projected images
n_images = 100; % number of images in simulation
n_r_perc = 50; % the radial resolution used in nufft as a percentage of image's size
max_shift_perc = 10; % the maximum shift to simulate in each image (x and y direction) as a percentage of image's size
shift_step = 0.5; % the shift step
mask_radius_perc = 50; % the mask applied to each image as a percentage of image's size
inplane_rot_res = 1;  % the resolution in angles for estimating the in-plane rotation angles
verbose = 1; % currently either 1 (detailed debug output), or 0 (minimal output debug messages)
output_folder = './results'; % folder where original and reconstructed volumes will reside

%% Step 0: Making sure that output folder exists
if exist(output_folder,'file') ~= 7
    error('Output folder %s does not exist. Please create it first.\n', output_folder);
end

%% Step 1: Creating cache file of all candidate rotation matrices in case n_symm > 4
% cache_file_name = 'cn_cache_points1000_ntheta360_res1.mat'; 
if n_symm > 4
    n_Points_sphere = 1000; % the number of viewing direction to sample (it is not recommended to use a number < 1000)
    log_message('Creating cache file under folder: %s',output_folder);
    log_message('#points on sphere=%d, n_theta=%d, inplane_rot_res=%d',n_Points_sphere,n_theta,inplane_rot_res);
    cache_file_name  = cryo_cn_create_cache(output_folder,n_Points_sphere,n_theta,inplane_rot_res);
end

%% Step 2: Simulate data %%
% A duplicated version (so that it is Cn-symmetric) of a non-symetric gaussian volume
max_shift = ceil(proj_size*max_shift_perc/100);
[projs,refq,~,~,vol_orig] = generate_cn_images(n_symm,n_images,snr,proj_size,'C1_DUPLICATED',max_shift,shift_step);
% saving original volume to disk
vol_orig_file_name = fullfile(output_folder,sprintf('vol_orig_c%d.mrc',n_symm));
log_message('saving original volume under %s',vol_orig_file_name);
WriteMRC(vol_orig,1,vol_orig_file_name);

if n_symm > 4
    % currently we don't know how to handle equators for n_symm > 4
    [projs,refq] = remove_eq_images(projs,refq);
end

if snr <= 1 % otherwise, masking makes thing worse
    mask_radius = proj_size*mask_radius_perc/100;
    log_message('Masking images using mask-radius=%d',mask_radius);
    masked_projs = mask_fuzzy(projs,mask_radius);
else
    masked_projs = projs;
    log_message('SNR=%.2e is greater than 1. Not performing mask', snr);
end

n_r = ceil(proj_size*n_r_perc/100);
log_message('Computing the polar Fourier transform of projections with resolution n_r=%d',n_r);
[npf,~] = cryo_pft(masked_projs,n_r,n_theta,'single');

if snr <= 1 % otherwise, filtering makes thing worse
    log_message('Guass filtering the images');
    npf = gaussian_filter_imgs(npf);
else
    log_message('SNR=%.2e is greater than 1. Not performing gauss filtering', snr);
end

%% Step 3: Computing the relative viewing directions
if(n_symm==3 || n_symm==4)
    max_shift_1d  = ceil(2*sqrt(2)*max_shift); % TODO: is this needed? if so, where?
    is_remove_non_rank1 = true; % whether or not to remove images with the least number of induced rank-1 matrices of relative viewing directions 
    non_rank1_remov_percent = 0.25; % the percent of images to remove
    log_message('Computing all relative viewing directions for n=3,4');
    [vijs,viis,npf,projs,refq] = compute_third_row_outer_prod_c34(n_symm,npf,max_shift_1d,shift_step,'',...
        projs,verbose,is_remove_non_rank1,non_rank1_remov_percent,refq);
else
    log_message('Computing all relative viewing directions for n>4');
    [vijs,viis,~] = compute_third_row_outer_prod_cn(npf,n_symm,max_shift,shift_step,cache_file_name,verbose,refq);
end

%% Step 4: Handedness synchronization
log_message('Handedness synchronization');
[vijs,viis] = global_sync_J(vijs,viis);

%% Step 5: Viewing direction estimation (i.e., estimating third row of each rotation matrix)
vis = estimate_third_rows(vijs,viis,n_symm);

%% step 6: In-plane rotations angles estimation
log_message('In plane angles estimation');
rots =      estimate_inplane_rotations(npf,vis,n_symm,inplane_rot_res,max_shift,shift_step);

%% Step 7: Results Analysis
[rot_alligned,err_in_degrees,mse] = analyze_results_cn(rots,n_symm,n_theta,refq);

%% Step 8: Reconstructing volume
log_message('Reconstructing volume');
estimatedVol = reconstruct_cn(projs,rot_alligned,n_symm,n_r,n_theta); % supply max_shift_2d or max_shift_1d ?

%% Step 9: Saving volumes
mrc_fname = sprintf('out_c%dnims%dshift%d.mrc',n_symm,n_images,max_shift_perc); 
recon_mrc_fname = fullfile(output_folder,mrc_fname); % name of reconstructed mrc volume
[vol_filename,vol_filt_file_name,vol_no_filt_file_name] = save_vols(estimatedVol,recon_mrc_fname,n_symm);