% Example code for generating class averages.
%
% The example generates 10000 noisy projection and produces 10000 class
% averages by averaging each image with its 50 nearest neighbors. 
%
% Data set:
% ---------
% 10000 noisy projections with additive white Gaussian noise at SNR/1=100.
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
% Tejal April 17, 2016
% New align_main: Tejal March 2017
fname='clean_data.mat';
if exist(fname,'file')==2
	data=load('clean_data.mat');
else
	gen_simulation_data;
end


K=10000; %K is the number of images
SNR = 1/10; %SNR
use_shifted=0;

use_CTF=1;
ndef=20; % Number of defocus groups
def1=1;
def2=4;
lambda = EWavelength(300);
B=10; % decay envelope parameter

[g_proj_CTF,CTF,defocus_group]=  add_CTF_env_v6(cfft2(data.projections(:,:,1:K)), ndef, def1,def2,B, lambda, use_CTF);
[images, noise_v_r]=addnoise_v6(icfft2(g_proj_CTF), SNR);
q = data.q(:, 1:K);
L = size(images, 1);
n_nbor = 10; %number of nearest neighbors for initial classification.
n_nbor_large=50;
isrann = 0;

use_VDM=1;
k_VDM_in = n_nbor; % number of nearest neighbors for building graph for VDM.
VDM_flag = 0; % VDM using union rule
k_VDM_out = n_nbor; % number of nearest neighbors search for 


% Initial classification with sPCA (new, fast code)
[ images_fl ] = Phase_Flip(images, defocus_group, CTF); %phase flipping 

str=which('gen_simulation_data.m');
[pathstr,~,~]=fileparts(str);
fname=fullfile(pathstr,'simulation','noisy_data_fl.mrc');
WriteMRC(real(images_fl),1,fname);
allims=imagestackReader(fname);

disp('Phase flipped');
[sPCA_data, sPCA_coeff_cell, basis, recon_spca]=data_sPCA(images_fl,  noise_v_r);
[mse_spca] = calc_MSE_v6(recon_spca, data.projections(:,:,1:K),sPCA_data.R);
tic_init=tic;
[ class_f, class_refl_f, rot_f, corr_f,  timing_f ] = Initial_classification_FD(sPCA_data, n_nbor, isrann );
toc_init=toc(tic_init);
disp('Finished initial classification...');
if(use_VDM)
	tic_VDM = tic;
	[ class_VDM, class_VDM_refl, rot_f_vdm ] = VDM(class_f, ones(size(class_f)), rot_f, class_refl_f, k_VDM_in, VDM_flag, k_VDM_out);
	toc_VDM = toc(tic_VDM);
	disp('Finished VDM classification...');
	% Check Classification result
	[ d_f, error_rot_f ] = check_simulation_results(class_VDM, class_VDM_refl, rot_f_vdm, q); % should use minus sign for init class, no minus sign for VDM 
else
	rot_f_vdm = rot_f;
	class_VDM = class_f;
	class_VDM_refl = class_refl_f;
	k_VDM_out = k_VDM_in;
	[ d_f, error_rot_f ] = check_simulation_results(class_f, class_refl_f, -rot_f, q); % should use minus sign for init class, no minus sign for VDM 
end

list_recon = [1:size(images_fl, 3)];

tmp_dir=tempmrcdir;
max_shift=0;
tic_align = tic;
[ shifts, corr, averagesfname, norm_variance ] = align_main( allims, rot_f_vdm, class_VDM, class_VDM_refl, sPCA_data, k_VDM_out, max_shift, list_recon, recon_spca, tmp_dir);
toc_align = toc(tic_align);
disp('Finished alignment and class averaging...');
% Check Classification result
d = d_f;
error_rot = error_rot_f;
[ N, X ] = hist(acosd(d), [0:180]);
figure; bar(N);
xlabel('a$\cos\langle v_i, v_j \rangle$', 'interpreter', 'latex');



