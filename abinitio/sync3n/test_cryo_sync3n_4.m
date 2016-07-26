% Test the 3Nx3N synchronization algorithm on noisy simulated data.
%
% Same as test_cryo_sync3n_3, but uses wieghts for the blocks of the
% synchronization matrix. See improvement in accuracy of the rotations over
% test_cryo_sync3n_3.
%
% Yoel Shkolnisky, April 2016.

K=400;
initstate;
q=qrand(K);

%L=1.0E15; % L is very large to reduce discretization errors due to finite angular resolution.
L=360;
clmatrix=clmatrix_cheat_q(q,L);
p=0.5;
%p=1;
[noisy_cl,is_perturbed]=perturb_clmatrix(clmatrix,L,p);

open_log(0);
use_weights=1;
rotations=cryo_sync3n_estimate_rotations(noisy_cl,L,use_weights);
dir1=Q2S2(q,L);
dir2=R2S2(rotations,L);
check_orientations(dir1,dir2);