tic;

initstate; 
open_log(0);

%% Load and display projections
% The MAT file p100_c4_shifted contains 100 projections of size 65x65. The
% orientations (given as quaternions) used to generate these projections
% are stored in the the variable "refq". The projection were generated using the following command:

max_shift  = 5; % number of shifts to consider
shift_step = 1; % the shift step (see Yoel's technical report)
[projs,refq] = generate_c2_images(100,10000000000,65,'GAUSSIAN',max_shift,shift_step);

% load p100_c2_gaussian_no_shifts;
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
min_dist_cls = 25; % the minimal distance (in degrees) between two lines in a single images
clmatrix = cryo_clmatrix_c2_gpu_tmp(npf,size(npf,3),1,max_shift,shift_step,min_dist_cls); 
cl_detection_rate(clmatrix,n_theta,refq);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 5  : calculate relative-rotations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Rijs,Rijgs,Confijs] = cryo_generateRij(clmatrix,n_theta,refq);
% Rijs  = cryo_c2_estimate_all_Rijs(clmatrix(:,:,1),n_theta,refq);
% Rijgs = cryo_c2_estimate_all_Rijs(clmatrix(:,:,2),n_theta,refq);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 6  : inner J-synchronization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nImages = size(clmatrix,1);
[Rijs,Rijgs] = local_sync_J(Rijs,Rijgs,nImages);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 7  : outer J-synchronization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Rijs,Rijgs] = global_sync_J(Rijs,Rijgs,nImages);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 8  : third rows estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vis  = estimate_third_rows(Rijs,Rijgs,nImages);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 9  : in-plane rotations angles estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rots = cryo_inplane_rotations(vis,Rijs,Rijgs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 10  : Results Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[rot_alligned,err_in_degrees,mse] = analyze_results(rots,n_theta,refq);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 11  : Reconstructing volume
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
estimatedVol = reconstruct(projs,rot_alligned,n_r,n_theta,max_shift,shift_step);   

% WriteMRC(estimatedVol,1,'example1_with_shifts.mrc');