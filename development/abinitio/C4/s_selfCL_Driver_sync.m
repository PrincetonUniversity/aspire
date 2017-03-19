tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

is_use_gpu = true;

max_shift  = 0; % number of shifts to consider
shift_step = 0.5; % the shift step (see Yoel's technical report)

n_theta = 360; % number of rays in every projection
n_r     = 89;  % number of radial points in every radial line

is_detect_equator_ims = true;

is_remove_non_rank1 = true;
non_rank1_remov_percent = 0.25; 

% the resolution in degrees of in-plane rotation angle. 
% resolution=1 means we look for angles 0,1,2,3 ...
inplane_rot_res = 0.25; 
is_reconstructVol = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           CODE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
initstate; 
open_log(0);
is_simulation = true;

if ~is_simulation
    noisy_projs = ReadMRC('/mnt/ix2/backup/datasets/10005/output/averages_nn100_group1.mrc');
else
    params_simul.c4_type = 'GAUSSIAN'; % either GAUSSIAN,SYNTHETIC or IP3 or TRPV1
    params_simul.nImages = 100;
    params_simul.SNR = 100000000;
%     params_simul.SNR = 1/16;
    params_simul.projSize = 65;  % Size of the projections. Can be even or odd.
    [noisy_projs,refq,ref_shifts] = generate_c4_images(params_simul.nImages,...
        params_simul.SNR,params_simul.projSize,params_simul.c4_type,max_shift,shift_step);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 1  : Computing polar Fourier transform of projections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~is_simulation
    error('need to specify radius for circular mask');
%     masked_projs = mask_fuzzy(noisy_projs,28);
else
    masked_projs = circular_mask(noisy_projs,params_simul.c4_type);
end

[npf,~] = cryo_pft(masked_projs,n_r,n_theta,'single');

npf = gaussian_filter_imgs(npf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 2  : detect a single pair of common-lines between each pair of images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,n_theta,nImages] = size(npf);
if is_use_gpu
    %TODO: why only integer shift_step? I think i need 0.1
    clmatrix = cryo_clmatrix_gpu(npf,nImages,1,max_shift,shift_step); 
else
    clmatrix = cryo_clmatrix(npf,nImages,1,max_shift,shift_step);
end

if is_simulation
    cl_detection_rate(clmatrix,n_theta,refq);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 3  : detect self-common-lines in each image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
equator_res_fact = 10;
equator_removal_frac = 0.1;
sclmatrix = cryo_self_clmatrix_gpu(npf,max_shift,shift_step,is_detect_equator_ims,equator_res_fact,...
    equator_removal_frac,refq);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 4  : calculate self-relative-rotations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Riis = estimate_all_Riis(sclmatrix,n_theta,refq);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 5  : calculate relative-rotations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rijs = cryo_c4_estimate_all_Rijs(clmatrix,n_theta,refq);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 6  : inner J-synchronization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[vijs,viis,npf,noisy_projs,refq,ref_shifts] = local_sync_J(Rijs,Riis,npf,...
                                noisy_projs,is_remove_non_rank1,non_rank1_remov_percent,refq,ref_shifts);

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
[rots,in_plane_rotations] = estimate_inplane_rotations2(npf,vis,inplane_rot_res,max_shift,shift_step);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 10  : Results Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if is_simulation
    [rot_alligned,err_in_degrees,mse] = analyze_results(rots,n_theta,refq);
else
    rot_alligned = rots;
    mse = Inf;
    err_in_degrees = Inf;
end
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 11  : Reconstructing volume
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if is_reconstructVol
    estimatedVol = reconstruct_vol(noisy_projs,npf,rot_alligned,max_shift,shift_step);
else
    estimatedVol = [];
end

run_time = toc