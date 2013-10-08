% Example code:
%
%   Example code for generating class averages.
%
% The example generates 10000 noisy projection and produces 10000 class
% averages by averaging each image with its 50 nearest neighbors. 
%
% Data set:
% ---------
% 10000 noisy projections of the 70S ribosomal subunit with additive white
% Gaussian at SNR/1=50. The 
% projections orientations are random and uniformly disributed, and they
% are shifted randomly with maximum shift of +/- 4 pixels in the x and y
% direcrtions. The clean projections are separated in 20 different defocus
% groups.
% If the file 'clean_data.mat' (in the folder ./simulation) is missing,
% please run  
%   gen_simulation_data 
% to generate it, before running this example code. 

clc;
clear all;
clf;
K = 10000; %K is the number of images
SNR = 1/50; %SNR
data = load('clean_data.mat'); % load clean centered projection images 
[images, defocus_group, noise, noise_spec, c, q]=create_image_wCTF_wshifts(data, SNR, 'gaussian'); %create projection images with CTF and shifts
clear data;
[ images ] = Phase_Flip(images, defocus_group, c); %phase flipping 
%Low pass filtering images to make it approximately invaraint to small
%shifts
[ images_lpf ] = low_pass_filter( images );
L = size(images, 1);
r_max = floor(L/2)-10;
n_nbor = 50;
isrann = 0;
k_VDM_in = 5; % number of nearest neighbors for building graph for VDM.
VDM_flag = 0;
k_VDM_out = 50; % number of nearest neighbors search for 
max_shift = 12;
% Initial Classification
[ class, class_refl, rot, ~, timing ] = Initial_classification(images_lpf, r_max, n_nbor, isrann );
clear images_lpf;
disp('Finished initial classification...');
% VDM Classification. Search for 50 nearest neighbors
[ class_VDM, class_VDM_refl, angle ] = VDM(class, ones(size(class)), rot, class_refl, k_VDM_in, VDM_flag, k_VDM_out);
disp('Finished VDM classification...');
% align and generate class averages
[ shifts, corr, average, norm_variance ] = align_main( images, angle, class_VDM, class_VDM_refl, k_VDM_out, max_shift);
disp('Finished alignment and class averaging...');
% Check Classification result
[ d, error_rot ] = check_simulation_results(class_VDM, class_VDM_refl, angle, q);
[ N, X ] = hist(acosd(d), [0:180]);
figure; bar(N);
xlabel('a$\cos\langle v_i, v_j \rangle$', 'interpreter', 'latex');
