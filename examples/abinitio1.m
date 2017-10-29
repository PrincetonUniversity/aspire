%% Basic abinitio reconstruction - example 1.
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

load p100_shifted;
viewstack(projections,10,10);   % Display the proejctions.

%% Compute reference common lines matrix
% Compute the true common lines matrix for the projections, by using their
% true orientations.

n_theta=72; % Angular resolution - number of sinograms computed for each 
            % projection. This corresponds to a resolution of 5 degrees.
[ref_clmatrix,~]=clmatrix_cheat(rots,n_theta);

%% Compute common lines from projections

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
[ clstack,~, shift_equations,~] = ...
    commonlines_gaussian(npf,max_shift,shift_step );

% Compare common lines computed from projections to the reference common
% lines. A common line is considered correctly-identified if it deviates
% from the true common line between the projections by up to 10 degrees.
prop=comparecl( clstack, ref_clmatrix, n_theta, 10 );
fprintf('Percentage of correct common lines: %f%%\n\n',prop*100);
[ est_shifts, shift_err] = test_shifts( shift_equations, shifts);
est_shifts=full(est_shifts);

%% Assign orientation using common lines, using least squares method.
% The resulting MSE should be small (of the order of 1e-4).

[est_inv_rots] = est_orientations_LUD(clstack, n_theta);
fprintf('MSE of the estimated rotations: %f\n\n', ...
    check_MSE(est_inv_rots,rots));


%% 3D inversion
[ v, v_b, kernel ,err, iter, flag] = recon3d_firm( projections,...
est_inv_rots,-est_shifts, 1e-6, 30, zeros(65,65,65));

assert(norm(imag(v(:)))/norm(v(:))<1.0e-4);
v=real(v);
WriteMRC(v,1,'example1.mrc'); % Output density map reconstructed from projections.
fprintf('Saved example1.mrc\n');

