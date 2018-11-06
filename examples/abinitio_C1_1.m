%% Basic abinitio reconstruction.
%
% The script demonstrates calling various functions of the aspire package.
%
% Follow comments below.
%
% Yoel Shkolnisky and Lanhui Wang, September 2013.

clear;

%% Load and display projections
% The MAT file p100_shifted contains 100 projections of size 65x65. The
% orientations (given as rotation matrices) used to generate these projections
% are stored in the the variable "rots". The shift introduced into each
% projection is stored in the variable "shifts". The projections were
% shifted by random integer shifts of up to +/- 5 pixels, with steps of 1
% pixel.

log_message('Loading and displaying images')
load p100_shifted;
viewstack(projections,10,10);   % Display the proejctions.

%% Compute reference common lines matrix
% Compute the true common lines matrix for the projections, by using their
% true orientations.

log_message('Computing ground-truth common lines');
n_theta=72; % Angular resolution - number of sinograms computed for each 
            % projection. This corresponds to a resolution of 5 degrees.
[ref_clmatrix,~]=clmatrix_cheat(rots,n_theta);

%% Compute common lines from projections

log_message('Estimating common lines from projections');
% Mask projections
mask_radius = 55;
[np,~]=mask_fuzzy(projections,mask_radius);

% Compute polar Fourier transform, using radial resolution n_r and angular
% resolution n_theta. n_theta is the same as above.
n_r=33;
[npf,~]=cryo_pft(np,n_r,n_theta);

% Find common lines between projections. Note that the max_shift used is
% much larger than above. This is because it is a 1D shift used to seearch
% for common lines. If d is the maximal expected 2D shift in a projection
% (in any direction), then the 1D shift used to find common lines should
% equal ceil(2*sqrt(2)*d).  
max_shift = 15;
shift_step = 1;
[clstack,~,shift_equations]= cryo_clmatrix(npf,-1,1,max_shift,shift_step);

% Compare common lines computed from projections to the reference common
% lines. A common line is considered correctly-identified if it deviates
% from the true common line between the projections by up to 10 degrees.
prop=comparecl( clstack, ref_clmatrix, n_theta, 10 );
fprintf('Percentage of correct common lines: %f%%\n\n',prop*100);
[ est_shifts, shift_err] = test_shifts( shift_equations, shifts);

%% Assign orientation using common lines, using least squares method.
% The resulting MSE should be small (of the order of 1e-4).

log_message('Estimating orientations of projections');
S=cryo_syncmatrix_vote(clstack,n_theta);
[est_inv_rots,diff,mse]=cryo_syncrotations(S,rots);
fprintf('MSE of the estimated rotations: %5.2e\n\n', ...
    check_MSE(est_inv_rots,rots));


%% 3D inversion
log_message('Reconstructing 3D density from projections');
params = struct();
params.rot_matrices = est_inv_rots;
params.ctf = ones(size(projections, 1)*ones(1, 2));
params.ctf_idx = ones(size(projections, 3), 1);
params.shifts = full(est_shifts');
params.ampl = ones(size(projections, 3), 1);

basis = dirac_basis(size(projections, 1)*ones(1, 3));

v = cryo_estimate_mean(projections, params, basis);

fname='example1.mrc';
WriteMRC(v,1,fname); % Output density map reconstructed from projections.
log_message('Reconstructed density saved to %s',fname);
log_message('Done!');
 
% load('cleanrib.mat','volref');
% WriteMRC(volref,1,'volref.mrc');
% cryo_compare_volumes('example1.mrc','volref.mrc',0.5,1,1);
% delete('volref.mrc');
