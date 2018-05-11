% Test the function cryo_sync3n_estimate_Rij_scores using noisy simulated
% data. 
%
% The script generates K quternions (images), constructs their
% corresponding common lines matrix, and introduces errors into the
% common lines matrix, such that only a factrion p of the common lines
% remain correct. 
% Then, the function calls cryo_sync3n_estimate_Rij_scores to estimate the
% percentage of correct common lines. The printed number P should
% comparable with p.
%
% Yoel Shkolnisky, March 2016.

K=400;
initstate;
rots = rand_rots(K);

%L=1.0E15; % L is very large to reduce discretization errors due to finite angular resolution.
L=360;
clmatrix=clmatrix_cheat(rots,L);
% XXX p=0.25 crashes!!!
p=0.3; % The estimation is accurate up p~0.2
%p=1;
[noisy_cl,is_perturbed]=perturb_clmatrix(clmatrix,L,p);

open_log(0);
[Rijs, ~, ~, ~] = cryo_sync3n_estimate_all_Rijs(noisy_cl, L);

required_P=0;
PLOT=1;
verbose=1;
[P, sigma, Rsquare, Pij, scores_hist, cum_scores, plotted] =...
    cryo_sync3n_estimate_Rij_scores(Rijs,required_P, PLOT, verbose);

log_message('True percentage of correct common lines: %5.2f',p);
log_message('Estimated percentage of correct common lines: %5.2f',P);
