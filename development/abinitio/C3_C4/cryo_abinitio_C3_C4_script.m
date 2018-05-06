% instack = '/scratch/yoel/denoised/10081/denoised_group1.mrcs';

n_symm = 3;
if n_symm == 3
    instack = '10004_averages_nn100_selected_small.mrcs';
    recon_folder = './results_tmp_C3';
else
    instack = '10081_denoised_group1_small.mrcs';
    recon_folder = './results_tmp_C4';
end

outvol = fullfile(recon_folder,'out.mrc');

%% option 1: using minimal set of input variables
cryo_abinitio_C3_C4(n_symm,instack,outvol);

%% option 2: using all input variables
outmat = fullfile(recon_folder,'out.mat');
max_shift_perc = 15;
shift_step = 0.5;
n_r_perc = 50;
mask_radius_perc = 70;
n_theta = 360;

log_fname = fullfile(recon_folder,'log.txt');
open_log(log_fname);
cryo_abinitio_C3_C4(n_symm,instack,outvol,outmat,max_shift_perc,shift_step,n_r_perc,mask_radius_perc,n_theta);
close_log();