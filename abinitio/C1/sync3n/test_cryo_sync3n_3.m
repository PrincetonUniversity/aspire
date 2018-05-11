% Test the 3Nx3N synchronization algorithm on noisy  simulated data.
%
% The script generates K rotations (images), constructs their
% corresponding common lines matrix, and introduces errors into the
% common lines matrix, such that only a factrion p of the common lines
% remain correct. Then, it estimates the rotations from the noisy matrix
% and compares the results to the true orientations. To display the error,
% the functions plots the histogram of the angle between the direction of
% each of the KL true Fourier lines and its estimated orientation.
%
% This functions is the same as ../sync2n/test_cryo_syncmatrix_3.m but uses
% the 3Nx3N algorithm instead of the 2Nx2N.
%
% Yoel Shkolnisky, March 2016.

K=400;
initstate;
rots_ref = rand_rots(K);

%L=1.0E15; % L is very large to reduce discretization errors due to finite angular resolution.
L=360;
clmatrix=clmatrix_cheat(rots_ref,L);
p=0.5;
%p=1;
[noisy_cl,is_perturbed]=perturb_clmatrix(clmatrix,L,p);

open_log(0);
rotations=cryo_sync3n_estimate_rotations(noisy_cl,L);
inv_rots_ref = permute(rots_ref, [2 1 3]);
dir1=R2S2(inv_rots_ref,L);
dir2=R2S2(rotations,L);
check_orientations(dir1,dir2);
