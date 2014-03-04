%test code
% If you don't have clean_data.mat in the folder ./simulation, please run
% gen_simulation_data before running this code.
% gen_simulation_data generates 10000 clean centered projection images.
%
%10^3 clean centered projection images.
% Variables determined by users are the following:
% r_max =floor(L/2); %radius of region of interest that contains the
% particle.
% n_nbor = 50; %number of nearest neighbors for initial classification.
% k_VDM_in = 5; % number of nearest neighbors for building graph for VDM.
% VDM_flag = 0; % Using union rule (0) or intersection rule (1) for
% constructring VDM matrix.
% k_VDM_out = 5; % output number of nearest neighbors
% max_shift = 12; % shift search range.

data = load('clean_data.mat');
K=1000;
images = data.projections(:, :, 1:K);
q = data.q(:, 1:K);

L = size(images, 1);
r_max = floor(L/2);
n_nbor = 50;
isrann = 0;
k_VDM_in = 5; % number of nearest neighbors for building graph for VDM.
VDM_flag = 0; % VDM using union rule
k_VDM_out = 5; % number of nearest neighbors search for 
max_shift = 12;
list_recon = 1:K;
%%low pass filtering the images (This step is optional for centered images. This makes the images roughly invariant to very small shifts)
%[ images_lpr ] = low_pass_filter(images);
% Initial Classification
[ class, class_refl, rot, ~, FBsPCA_data, timing ] = Initial_classification(images, r_max, n_nbor, isrann );
disp('Finished initial classification...');
% VDM Classification. Search for 50 nearest neighbors
tic_VDM = tic;
[ class_VDM, class_VDM_refl, angle ] = VDM(class, ones(size(class)), rot, class_refl, k_VDM_in, VDM_flag, k_VDM_out);
toc_VDM = toc(tic_VDM);
disp('Finished VDM classification...');
% align and generate class averages
tic_align = tic;
[ shifts, corr, average, norm_variance ] = align_main( images, angle, class_VDM, class_VDM_refl, FBsPCA_data, k_VDM_out, max_shift, list_recon );
toc_align = toc(tic_align);
disp('Finished alignment and class averaging...');
% Check Classification result
[ d, error_rot ] = check_simulation_results(class_VDM, class_VDM_refl, angle, q);
[ N, X ] = hist(acosd(d), [0:180]);
figure; bar(N);
xlabel('a$\cos\langle v_i, v_j \rangle$', 'interpreter', 'Latex');