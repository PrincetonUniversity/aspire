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
output_folder = fullfile(pwd,'results'); % folder where original and reconstructed volumes will reside

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
else
    cache_file_name = ''; % cache is only needed in case n>4
end

%% Step 2: Simulate data %%
% A duplicated version (so that it is Cn-symmetric) of a non-symetric gaussian volume
max_shift = ceil(proj_size*max_shift_perc/100);
[projs,refq,~,~,vol_orig] = generate_cn_images(n_symm,n_images,snr,proj_size,'C1_DUPLICATED',max_shift,shift_step);

if n_symm > 4
    % currently we don't know how to handle equators for n_symm > 4
    [projs,refq] = remove_eq_images(projs,refq);
end

instack = fullfile(output_folder,sprintf('simulated_instack_c%d.mrc',n_symm));
log_message('Saving simulated projection images under %s',instack);
WriteMRC(projs,1,instack);

% saving original volume to disk so that we can later compare the results
vol_orig_mrc_name = fullfile(output_folder,sprintf('vol_orig_c%d.mrc',n_symm));
log_message('Saving original volume under %s',vol_orig_mrc_name);
WriteMRC(vol_orig,1,vol_orig_mrc_name);

outvol = sprintf('out_c%dnims%dshift%d.mrc',n_symm,n_images,max_shift_perc); 
outvol = fullfile(output_folder,outvol); % name of reconstructed mrc volume

matfname = ''; % supply mat file name if you wish to save intermediate results 


%% Step 3: Call main impl function %%
cryo_abinitio_Cn(n_symm,instack,outvol,cache_file_name,matfname,...
    verbose,n_theta,n_r_perc,max_shift_perc,shift_step,mask_radius_perc,inplane_rot_res,refq);

% if snr <= 1 % otherwise, masking makes thing worse
%     mask_radius = proj_size*mask_radius_perc/100;
%     log_message('Masking images using mask-radius=%d',mask_radius);
%     masked_projs = mask_fuzzy(projs,mask_radius);
% else
%     masked_projs = projs;
%     log_message('SNR=%.2e is greater than 1. Not performing mask', snr);
% end
% 
% n_r = ceil(proj_size*n_r_perc/100);
% log_message('Computing the polar Fourier transform of projections with resolution n_r=%d',n_r);
% [npf,~] = cryo_pft(masked_projs,n_r,n_theta,'single');
% 
% if snr <= 1 % otherwise, filtering makes thing worse
%     log_message('Guass filtering the images');
%     npf = gaussian_filter_imgs(npf);
% else
%     log_message('SNR=%.2e is greater than 1. Not performing gauss filtering', snr);
% end
% 