tic;

initstate; 
open_log(0);

%% Load and display projections
% The MAT file p100_c4_shifted contains 100 projections of size 65x65. The
% orientations (given as quaternions) used to generate these projections
% are stored in the the variable "refq". The shift introduced into each
% projection is stored in the variable "ref_shifts". The projections were
% shifted by random integer shifts of up to +/- 5 pixels, with steps of 0.5
% pixel.
load p100_c4_shifted;
viewstack(projs,5,5);   % Display the proejctions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 1  : Computing polar Fourier transform of projections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
masked_projs = mask_fuzzy(projs,23);

n_theta = 360; % number of rays in every projection
n_r     = 89;  % number of radial points in every radial line
[npf,~] = cryo_pft(masked_projs,n_r,n_theta,'single');

npf = gaussian_filter_imgs(npf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 2  : detect a single pair of common-lines between each pair of images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
max_shift  = 5; % number of shifts to consider
shift_step = 0.5; % the shift step (see Yoel's technical report)
clmatrix = cryo_clmatrix_gpu(npf,size(npf,3),1,max_shift,shift_step); 
cl_detection_rate(clmatrix,n_theta,refq);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 3  : detect self-common-lines in each image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
is_handle_equator_ims = true;
equator_res_fact = 10;
equator_removal_frac = 0.1;
sclmatrix = cryo_self_clmatrix_gpu(npf,max_shift,shift_step,...
    is_handle_equator_ims,equator_res_fact,equator_removal_frac,refq);

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
is_remove_non_rank1 = true;
non_rank1_remov_percent = 0.25;
[vijs,viis,npf,projs,refq,ref_shifts] = local_sync_J(Rijs,Riis,npf,...
                                projs,is_remove_non_rank1,non_rank1_remov_percent,refq,ref_shifts);

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
inplane_rot_res = 1;
[rots,in_plane_rotations] = estimate_inplane_rotations2(npf,vis,inplane_rot_res,max_shift,shift_step);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 10  : Results Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[rot_alligned,err_in_degrees,mse] = analyze_results(rots,n_theta,refq);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 11  : Reconstructing volume
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
estimatedVol = reconstruct_vol(projs,npf,rot_alligned,max_shift,shift_step);