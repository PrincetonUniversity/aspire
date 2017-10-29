% Example code for generating class averages.
%
% The example generates 10000 noisy projection and produces 10000 class
% averages by averaging each image with its 50 nearest neighbors. 
%
% Data set:
% ---------
% 1000 noisy projections with additive white Gaussian noise at SNR/1=100.
% The 
% projections orientations are random and uniformly disributed, and they
% are shifted randomly with maximum shift of +/- 4 pixels in the x and y
% direcrtions. The clean projections are separated in 20 different defocus
% groups.
% If the file 'clean_data.mat' (in the folder ./simulation) is missing,
% please go into "ASPIRE/projections/class_average/" and run  
%   gen_simulation_data 
% to generate it, before running this example code. 
%
%
% Tuning the class averaging algorithm:
% -------------------------------------
% Variables determined by the user are the following
%
%   r_max =floor(L/2)-10; 
%       Radius of region of interest that contains the particle.
%
%   n_nbor = 50; 
%       Number of nearest neighbors for initial classification.
%
% Use new FB coefficients, new Init_class function
% Tejal April 17, 2016, modified Zhizhen Zhao Sep 2017: remove steerable PCA reconstruction
% New align_main: Tejal March 2017, Zhao align_main: modified September 2017
fname='clean_data.mat';
if exist(fname,'file')~=2
    log_message('Dataset not found. Generating...(may take time)');
    gen_simulation_data;
    log_message('Finished generating data');
end

log_message('Loading dataset %s',fname);
data=load('clean_data.mat');
log_message('Finished loading dataset');

K=1000; %K is the number of images
SNR = 1/10; %SNR
%SNR=1000000;
use_shifted=0;
log_message('Using K=%d images with use_shifts=%d',K,use_shifted);

use_CTF=1;
ndef=20; % Number of defocus groups
min_defocus=1;
max_defocus=4;
lambda = EWavelength(300);
B=10; % decay envelope parameter

log_message('Applying CTF to images')
log_message('Number of defocus groups is %d. min_defocus=%d, max defocus=%d',...
    ndef,min_defocus,max_defocus);
[g_proj_CTF,CTF,defocus_group]=  add_CTF_env_v6(cfft2(data.projections(:,:,1:K)), ndef, min_defocus,max_defocus,B, lambda, use_CTF);
log_message('Finished applying CTF');
log_message('Adding noise to images at SNR=%d',SNR)
[images, noise_v_r]=addnoise_v6(icfft2(g_proj_CTF), SNR);
log_message('Finished adding noise');

clear g_proj_CTF
rots_ref = data.rots(:, :, 1:K);
L = size(images, 1);
n_nbor = 10; %number of nearest neighbors for initial classification.
n_nbor_large=50;
isrann = 0;

use_VDM=1;
k_VDM_in = n_nbor; % number of nearest neighbors for building graph for VDM.
VDM_flag = 0; % VDM using union rule
k_VDM_out = n_nbor; % number of nearest neighbors search for 

% Initial classification with sPCA (new, fast code)
log_message('Phase flipping images')
images_fl = Phase_Flip(images, defocus_group, CTF); %phase flipping 
log_message('Finished phase flipping');

log_message('Save phase flipped images');
clear images 
str=which('gen_simulation_data.m');
[pathstr,~,~]=fileparts(str);
fname=fullfile(pathstr,'simulation','noisy_data_fl.mrc');
WriteMRC(real(images_fl),1,fname);
allims=imagestackReader(fname);
log_message('Finished saving');

log_message('Starting fast steerable PCA to compress and denoise images')
[sPCA_data]=data_sPCA(images_fl,  noise_v_r);
%[mse_spca] = calc_MSE_v6(recon_spca, data.projections(:,:,1:K),sPCA_data.R);
log_message('Finished steerable PCA');
%log_message('Relative MSE of denoised images after PCA is %f',mse_spca)

log_message('Starting initial classification');
tic_init=tic;
[ class_f, class_refl_f, rot_f, corr_f,  timing_f ] = Initial_classification_FD_update(sPCA_data, n_nbor, isrann );
toc_init=toc(tic_init);
log_message('Finished initial classification');
disp('Improving initial classification with vector diffusion maps...');

if(use_VDM)
    log_message('Starting VDM classification');
	tic_VDM = tic;
	[ class_VDM, class_VDM_refl, rot_f_vdm ] = VDM(class_f, ones(size(class_f)), rot_f, class_refl_f, k_VDM_in, VDM_flag, k_VDM_out);
	toc_VDM = toc(tic_VDM);
	log_message('Finished VDM classification');
	% Check Classification result
    log_message('Printing simuation results with VDM classification');
	[ d_f, error_rot_f ] = check_simulation_results(class_VDM, class_VDM_refl, rot_f_vdm, rots_ref); % should use minus sign for init class, no minus sign for VDM 
else
	rot_f_vdm = rot_f;
	class_VDM = class_f;
	class_VDM_refl = class_refl_f;
	k_VDM_out = k_VDM_in;
    log_message('Printing simuation results WITHOUT VDM classification');
	[ d_f, error_rot_f ] = check_simulation_results(class_f, class_refl_f, -rot_f, rots_ref); % should use minus sign for init class, no minus sign for VDM 
end

log_message('Start generating class averages (align_main)');
list_recon = 1:size(images_fl, 3);
%[tmp_dir]=fmtinput('Please enter the destination path to an empty folder for class averaged images. The default destination is your home directory','~/.','%s');
tmp_dir=tempmrcdir;
log_message('Using temporary folder %s',tmp_dir);
max_shift=0;
log_message('Starting align_main');
tic_align = tic;
[ shifts, corr, averagesfname, norm_variance ] = align_main( allims, rot_f_vdm, class_VDM, class_VDM_refl, sPCA_data, k_VDM_out, max_shift, list_recon, tmp_dir);
toc_align = toc(tic_align);
log_message('Finished generating class averages');
log_message('Checking simulation results');
% Check Classification result
d = d_f;
error_rot = error_rot_f;
[ N, X ] = hist(acosd(d), 0:180);
figure; bar(N);
xlabel('a$\cos\langle v_i, v_j \rangle$', 'interpreter', 'latex');
log_message('The histogram shows the angular distance in degrees between images classified into the same class. In the case of good classification, this distance should be as close to zero as possible.')
log_message('Errors (degrees): mean=%5.2f, median=%5.2f, std=%5.2f',...
    mean(acosd(d)),median(acosd(d)), std(acosd(d)));
log_message('Removing temporary files')
delete(sprintf('%s/*',tmp_dir));
rmdir(tmp_dir);
log_message('Test done!');
