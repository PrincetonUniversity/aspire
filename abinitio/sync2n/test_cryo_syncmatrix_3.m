% Test the synchronization algorithm for noisy  simulated data.
%
% The script generates K quternions (images), constructs their
% corresponding common lines matrix, and introduces errors into the
% common lines matrix, such that only a factrion p of the common lines
% remain correct. Then, it estimates the rotations from the noisy matrix
% and compares the results to the true orientations. To display the error,
% the functions plots the histogram of the angle between the direction of
% each of the KL true Fourier lines and its estimated orientation.
%
% Yoel Shkolnisky, August 2010.
% Revision: Renamed from test7. Y.S. September 2013.


K=400;
initstate;
q=qrand(K);

%L=1.0E15; % L is very large to reduce discretization errors due to finite angular resolution.
L=360;
clmatrix=clmatrix_cheat_q(q,L);
p=0.5;
%p=1;
[noisy_cl,is_perturbed]=perturb_clmatrix(clmatrix,L,p);

%S=cryo_syncmatrix_vote(noisy_cl,L,q,is_perturbed);
S=cryo_syncmatrix_vote(noisy_cl,L);
rotations=cryo_syncrotations(S,q);
dir1=Q2S2(q,L);
dir2=R2S2(rotations,L);
check_orientations(dir1,dir2)
