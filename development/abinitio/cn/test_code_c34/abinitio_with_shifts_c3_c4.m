tic;

initstate; 
open_log(0);

n_symm = 3;

if n_symm ~= 3 && n_symm ~= 4
    error('n_symm may be either 3 or 4');
end


%% Load and display projections
% The MAT file p100_c4_shifted contains 100 projections of size 65x65. The
% orientations (given as quaternions) used to generate these projections
% are stored in the the variable "refq". The projection were generated using the following command:
max_shift_2d  = 5;
shift_step_2d = 1;


[projs,refq,~,~,vol_orig] = generate_cn_images(n_symm,100,1000000,65,'C1_Eytan',0,1);
% 
% if n_symm == 3 
%     [projs,refq] = generate_c3_images(100,1000000,65,'GAUSSIAN',max_shift_2d,shift_step_2d);
% else
%     [projs,refq] = generate_c4_images(100,1000000,65,'GAUSSIAN',max_shift_2d,shift_step_2d);
% end

% load p100_c4_gaussian_no_shifts;
viewstack(projs,5,5);   % Display the proejctions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 1  : Computing polar Fourier transform of projections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
masked_projs = mask_fuzzy(projs,50);

n_theta = 360; % number of rays in every projection
n_r     = 89;  % number of radial points in every radial line
[npf,~] = cryo_pft(masked_projs,n_r,n_theta,'single');

% npf = gaussian_filter_imgs(npf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 2  : detect a single pair of common-lines between each pair of images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
max_shift_1d  = ceil(2*sqrt(2)*max_shift_2d);
shift_step_1d = 0.5;
clmatrix = cryo_clmatrix_gpu(npf,size(npf,3),1,max_shift_1d,shift_step_1d); 
cl_detection_rate_c3_c4(n_symm,clmatrix,n_theta,refq);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 3  : detect self-common-lines in each image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if n_symm == 3
    is_handle_equator_ims = false;
else
    is_handle_equator_ims = true;
end
sclmatrix = cryo_self_clmatrix_gpu_c3_c4(n_symm,npf,max_shift_1d,shift_step_1d,is_handle_equator_ims,refq);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 4  : calculate self-relative-rotations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Riis = estimate_all_Riis_c3_c4(n_symm,sclmatrix,n_theta,refq);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 5  : calculate relative-rotations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rijs = cryo_estimate_all_Rijs_c3_c4(n_symm,clmatrix,n_theta,refq);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 6  : inner J-synchronization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
is_remove_non_rank1 = true;
non_rank1_remov_percent = 0.25;
[vijs,viis,im_inds_to_remove,pairwise_inds_to_remove,...
    npf,projs,refq] = local_sync_J_c3_c4(n_symm,Rijs,Riis,npf,...
                                projs,is_remove_non_rank1,non_rank1_remov_percent,refq);
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
[rots,in_plane_rotations] =      estimate_inplane_rotations(npf,vis,n_symm,inplane_rot_res,max_shift_1d,shift_step_1d);
% [rots,in_plane_rotations] = estimate_inplane_rotations2_c3_c4(n_symm,npf,vis,inplane_rot_res,max_shift_1d,shift_step_1d);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 10  : Results Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [rot_alligned,err_in_degrees,mse] = analyze_results_c3_c4(n_symm,rots,n_theta,refq);
[rot_alligned,err_in_degrees,mse] = analyze_results_ml(rots,n_symm,n_theta,refq);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 11  : Reconstructing volume
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
estimatedVol = reconstruct_c3_c4(n_symm,projs,rot_alligned,n_r,n_theta);   

WriteMRC(estimatedVol,1,'example1_with_shifts.mrc');