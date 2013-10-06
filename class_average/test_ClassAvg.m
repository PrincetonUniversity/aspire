%test code
%Generate 1000 clean centered projection images.
initstate;
K = 1000; %K is the number of images
q = qrand(K);
load volume.mat
images=cryo_project(vol,q);
%save /simulation/clean_data_1000 images q

L = size(images, 1);
r_max = floor(L/2);
n_nbor = 10;
isrann = 0;
k_VDM_in = 5; % number of nearest neighbors for building graph for VDM.
VDM_flag = 0; % VDM using union rule
k_VDM_out = 5; % number of nearest neighbors search for 
max_shift = 12;
% Initial Classification
[ class, class_refl, rot, ~, timing ] = Initial_classification(images, r_max, n_nbor, isrann );
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
xlabel('a$\cos\langle v_i, v_j \rangle$', 'interpreter', 'Latex');
