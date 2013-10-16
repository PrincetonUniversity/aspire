%Before running the test code in refinement, please go to class_average
%folder and check if you have the mat file 'clean_data.mat' with 10^4 projection images
%in the sub-folder ./simulation. CTF information are in
%create_image_wCTF_wshifts.m
%
%The reference volume (70S ribosome) is in aspire/class_average/simulation/volume.mat. The
%reference volume is of 129x129x129 voxels. 

data = load('clean_data.mat'); % load clean centered projection images 
SNR = 1;
[images, defocus_group, noise, noise_spec, c, q]=create_image_wCTF_wshifts(data, SNR, 'gaussian'); %create projection images with CTF and shifts
load volume
L=size(data, 1); %Note image size L should be odd

%CTF Parameters
params.c = c ; %
params.N=20; %20 defocus groups
params.defidx=defocus_group;
params.max_shifts=6;  %maximum shift search range

filename = 'refined_model';
iter_max = 2;
tol = 0.1;

[ v_new ] = Refine( vol, images, params, 1, 0.1, filename );