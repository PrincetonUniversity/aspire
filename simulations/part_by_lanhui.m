K=100;
n_r=100;
n_theta=72;

% load projections
max_shift=5;
step_size=1;
noise_type='gaussian';
silent=1;
mask_radius = 55;
load p100_shifted;

[ref_clmatrix,clcorr]=clmatrix_cheat_q(q,n_theta);
% Mask
[np,sigma]=mask_fuzzy(projections,mask_radius);
% PFT
[npf,freqs]=cryo_pft(np,n_r,n_theta);
max_shift = 15;
shift_step = 1;
% Commonline detection
[ clstack,corrstack, shift_equations,shift_equations_map]...
    = commonlines_gaussian( npf,max_shift,shift_step );
prop=comparecl( clstack, ref_clmatrix, n_theta, 10 );
fprintf('correctness of common lines: %f\n\n',prop);
[ est_shifts, err] = test_shifts( shift_equations, shifts);
% Orientation determination
[est_inv_rots] = est_orientations_LS(clstack, n_theta);
fprintf('MSE of the estimated rotations: %f\n\n',check_MSE(est_inv_rots,q));
% 3D inversion
[ v, v_b, kernel ,err, iter, flag] = recon3d_firm( projections,...
est_inv_rots,-est_shifts, 1e-6, 30, zeros(129,129,129));