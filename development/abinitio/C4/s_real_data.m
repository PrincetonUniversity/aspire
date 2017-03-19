tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_theta = 360; % number of rays in every projection
n_r     = 89;  % number of radial points in every radial line
is_reconstructVol = true;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           CODE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
initstate; 
open_log(0);

projs = ReadMRC('/mnt/ix2/backup/datasets/10005/output/averages_nn100_group1.mrc');
projs(:,:,5001:end) = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 1  : Computing polar Fourier transform of projections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
masked_projs = mask_fuzzy(projs,28);

[npf,~] = cryo_pft(masked_projs,n_r,n_theta,'single');

npf = gaussian_filter_imgs(npf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 2  : detect a single pair of common-lines between each pair of images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,~,nImages] = size(npf);
log_message('detecting common-lines');
clmatrix = cryo_clmatrix_gpu(npf,nImages,1); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 3  : detect self-common-lines in each image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sclmatrix = cryo_self_clmatrix_gpu(npf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 4  : calculate self-relative-rotations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Riis = estimate_all_Riis(sclmatrix,n_theta);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 5  : calculate relative-rotations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rijs = cryo_c4_estimate_all_Rijs(clmatrix,n_theta);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 6  : inner J-synchronization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[vijs,viis,npf,projs] = local_sync_J(Rijs,Riis,npf,projs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 7  : outer J-synchronization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[vijs,viis] = global_sync_J(vijs,viis);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 8  : third rows estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vis  = estimate_third_rows(vijs,viis);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 9  : in-plane rotations angles estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[rots,in_plane_rotations] = estimate_inplane_rotations2(npf,vis);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 10  : Results Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rot_alligned = rots;
mse = Inf;
err_in_degrees = Inf;

%% TODO : add this code
% 
% pixA = 3.49;
% 
% vol1 = ReadMRC('ip3_nn_30_group1_randims.mrc');
% vol2 = ReadMRC('ip3_nn_30_group2_randims.mrc');
% 
% [Rest,estdx,vol2aligned] = cryo_align_densities(vol1,vol2,0,1,[],0,50);
% 
% plotFSC(vol1,vol2aligned,0.143,pixA);
% 
% 
% [h1,h2]=cryo_plot_viewing_directions(rot_alligned);
% [h1,h2]=cryo_plot_viewing_directions(rots2);
%

run_time = toc

if is_reconstructVol
    rot_alligned = rots;
    estimatedVol = reconstruct_vol(projs,npf,rot_alligned);
end