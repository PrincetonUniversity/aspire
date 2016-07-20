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

clear all;
K = 1000; %K is the number of images
SNR = 1/10; %SNR
data = load('clean_data.mat'); % load clean centered projection images 
disp('Loaded clean data')
downsampleddim=65;
sprintf('Downsampling to %dX%d grid', downsampleddim, downsampleddim)
data.projections=cryo_downsample(data.projections,[downsampleddim downsampleddim],1);
use_CTF=1;
ndef=20; % Number of defocus groups
def1=1;
def2=4;
lambda = EWavelength(300);
B=10; % decay envelope parameter

[g_proj_CTF,CTF,defocus_group]=  add_CTF_env_v6(cfft2(data.projections(:,:,1:K)), ndef, def1,def2,B, lambda, use_CTF);
[images, noise_v_r]=addnoise_v6(icfft2(g_proj_CTF), SNR);
% Use this for test with clean data
q = data.q(:, 1:K);
%images = data.projections(:, :, 1:K);
%clear data;
L = size(images, 1);
n_nbor = 30; %number of nearest neighbors for initial classification.
isrann = 0;

% Initial classification with sPCA (new, fast code)
[ images ] = Phase_Flip(images, defocus_group, CTF); %phase flipping 
disp('Phase flipped');
[sPCA_data, sPCA_coeff_cell, basis, recon_spca]=data_sPCA(images,  noise_v_r);
[mse_spca] = calc_MSE_v6(recon_spca, data.projections(:,:,1:K),sPCA_data.R);
[ class_f, class_refl_f, rot_f, corr_f,  timing_f ] = Initial_classification_FD(sPCA_data, n_nbor, isrann );
disp('Finished initial classification...');
% Check Classification result
[ d_f, error_rot_f ] = check_simulation_results(class_f, class_refl_f, -rot_f, q); % should use minus sign for init class, no minus sign for VDM 
[ N_f, X_f ] = hist(acosd(d_f), [0:180]);
figure; bar(N_f); title('sPCA')
xlabel('a$\cos\langle v_i, v_j \rangle$', 'interpreter', 'latex');
sprintf('sPCA: Number of images with correlation > %f is %d',0.9, numel(find(d_f(d_f>=0.9))))


