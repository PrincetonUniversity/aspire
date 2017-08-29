% Test the 3Nx3N synchronization algorithm on noisy simulated data.
%
% Same as test_cryo_sync3n_3, but uses wieghts for the blocks of the
% synchronization matrix. See improvement in accuracy of the rotations over
% test_cryo_sync3n_3.
%
% Yoel Shkolnisky, April 2016.

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
use_weights=1;
rotations=cryo_sync3n_estimate_rotations(noisy_cl,L,use_weights);
inv_rots_ref = permute(rots_ref, [2 1 3]);
dir1=R2S2(inv_rots_ref,L);
dir2=R2S2(rotations,L);
check_orientations(dir1,dir2);
