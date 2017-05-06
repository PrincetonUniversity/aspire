load vol_nn100_grp1_3000;

tic;
rots = estimate_inplane_rotations4(npf,vis,1,max_shift,shift_step);
toc

clear;
load vol_nn100_grp1_3000;

tic;
rots = estimate_inplane_rotations2(npf,vis,1,max_shift,shift_step);
toc
