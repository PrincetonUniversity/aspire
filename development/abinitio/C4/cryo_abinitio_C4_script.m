% instack = '/scratch/yoel/denoised/10081/denoised_group1.mrcs';
instack = '10081_denoised_group1_small.mrcs';
recon_folder = 'results_tmp';
outvol = fullfile(recon_folder,'out.mrc');

%% option 1: using minimal set of input variables
cryo_abinitio_C4(instack,outvol);

%% option 2: using all input variables
outmat = fullfile(recon_folder,'out.mat');
max_shift_perc = 15;
shift_step = 0.5;
n_r_perc = 50;
mask_radius_perc = 70;
n_theta = 360;

log_fname = fullfile(recon_folder,'log.txt');
open_log(log_fname);
cryo_abinitio_C4(instack,outvol,outmat,max_shift_perc,shift_step,n_r_perc,mask_radius_perc,n_theta);
close_log();