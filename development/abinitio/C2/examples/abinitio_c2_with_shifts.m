tic;

initstate; 
open_log(0);

%% Load and display projections
% The MAT file p100_c4_shifted contains 100 projections of size 65x65. The
% orientations (given as quaternions) used to generate these projections
% are stored in the the variable "refq". The projection were generated using the following command:

max_shift = 5; % number of shifts to consider
shift_step_2d = 1; % the shift step (see Yoel's technical report)
shift_step_1d = 0.5;
nImages = 100;
projSize = 89;

[projs,refq] = generate_c2_images(nImages,1/8,projSize,'SYNTHETIC',max_shift,shift_step_2d);

% load p100_c2_gaussian_no_shifts
viewstack(projs,5,5);   % Display the proejctions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 1  : Computing polar Fourier transform of projections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% masked_projs = mask_fuzzy(projs,23);

n_theta = 360; % number of rays in every projection
n_r     = ceil(projSize/2);  % number of radial points in every radial line
[npf,~] = cryo_pft(projs,n_r,n_theta,'single');

npf = gaussian_filter_imgs(npf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 2  : detect a single pair of common-lines between each pair of images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
min_dist_cls = 30; % the minimal distance (in degrees) between two lines in a single images
% the maximum shift occurs at a diagonal direction
max_1d_shift = ceil(2*sqrt(2)*max_shift);
[clmatrix,corrstack] = cryo_clmatrix_c2_gpu(npf,nImages,1,max_1d_shift,shift_step_1d,min_dist_cls);
cl_detection_rate(clmatrix,n_theta,min_dist_cls,refq);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 5  : calculate relative-rotations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Rijs,Rijgs,confijs] = cryo_generateRij(clmatrix,n_theta,refq);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 6  : inner J-synchronization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nImages = size(clmatrix,1);
[Rijs,Rijgs,~,isRank1_ijs] = local_sync_J(Rijs,Rijgs,nImages);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 7  : outer J-synchronization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Rijs,Rijgs] = global_sync_J(Rijs,Rijgs,nImages);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 8  : third rows estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

conf = confijs.*isRank1_ijs; 
is_use_weights = true;
vis  = estimate_third_rows(Rijs,Rijgs,conf,nImages,is_use_weights);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 9  : in-plane rotations angles estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
is_use_inplane_weights = true;
rots = cryo_inplane_rotations(vis,Rijs,Rijgs,is_use_inplane_weights,conf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 10  : Results Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[rot_alligned,err_in_degrees,mse] = analyze_results(rots,n_theta,refq);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 11  : Reconstructing volume
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% estimatedVol = reconstruct(projs,rot_alligned,n_r,n_theta,max_shift,shift_step_2d);   

% WriteMRC(estimatedVol,1,'example1_with_shifts.mrc');