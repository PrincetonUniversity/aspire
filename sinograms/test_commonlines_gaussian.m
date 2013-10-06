%% Search common lines for images with shifts using Gaussian Filter
% Test on clean images
fprintf('Commonline detection for clean shifted images using Gaussian filter.\n');
fprintf('=====================================================================\n');
K=100;
n_r=100;
n_theta=72;

max_shift=5;
step_size=1;
noise_type='gaussian';
silent=1;
mask_radius = 55;
% fprecomp='p100_shifted';
% [p, np, shifts, q] = ...
%     gen_projections(K,1,max_shift,step_size,noise_type,silent,fprecomp);
[p, np, shifts, q] = ...
    gen_projections(K,1,max_shift,step_size,noise_type,silent);
[ref_clmatrix,clcorr]=clmatrix_cheat_q(q,n_theta);
[np,sigma]=mask_fuzzy(np,mask_radius);
[npf,freqs]=cryo_pft(np,n_r,n_theta);
max_shift = 15;
shift_step = 1;
[ clstack,corrstack, shift_equations,shift_equations_map]...
    = commonlines_gaussian( npf,max_shift,shift_step );
prop=comparecl( clstack, ref_clmatrix, n_theta, 10 );
fprintf('correctness of common lines: %f\n\n',prop);
[ est_shifts, err] = test_shifts( shift_equations, shifts);
hist(err)
title('L2 norm errors of the estimated shifts');
    

% Test on noisy images
fprintf('Commonline detection for noisy shifted images using Gaussian filter.\n');
fprintf('=====================================================================\n');
K=100;
n_r=100;
n_theta=72;
for k=1:6
    SNR=1/2^k;
    fprintf('\nSNR=1/%d, ',2^k)
    max_shift=5;
    step_size=1;
    noise_type='gaussian';
    silent=1;
    mask_radius = 55;
%    fprecomp='p100_shifted';
    [p, np, shifts, q] = ...
        gen_projections(K,SNR,max_shift,step_size,noise_type,silent);
    [ref_clmatrix,clcorr]=clmatrix_cheat_q(q,n_theta);
    [np,sigma]=mask_fuzzy(np,mask_radius);
    [npf,freqs]=cryo_pft(np,n_r,n_theta);
    max_shift = 15;
    shift_step = 1;
    [ clstack,corrstack, shift_equations,shift_equations_map]...
                        = commonlines_gaussian( npf,max_shift,shift_step );
    prop=comparecl( clstack, ref_clmatrix, n_theta, 10 );
    fprintf('correctness of common lines: %f\n\n',prop);
    fprintf('-----------------------------------------------------------------\n');
end

% Results:
% Commonline detection for shifted images using Gaussian filter.
% =====================================================================
% 
% SNR=1/2, correctness of common lines: 0.944646
% 
% -----------------------------------------------------------------
% 
% SNR=1/4, correctness of common lines: 0.833939
% 
% -----------------------------------------------------------------
% 
% SNR=1/8, correctness of common lines: 0.625253
% 
% -----------------------------------------------------------------
% 
% SNR=1/16, correctness of common lines: 0.397374
% 
% -----------------------------------------------------------------
% 
% SNR=1/32, correctness of common lines: 0.212727
% 
% -----------------------------------------------------------------
% 
% SNR=1/64, correctness of common lines: 0.091515
% 
% -----------------------------------------------------------------
