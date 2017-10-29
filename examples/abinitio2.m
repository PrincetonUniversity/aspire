% Basic abinitio reconstruction - example 2.
%
% The script demonstrates calling various functions of the aspire package.
%
% Follow comments below.
%
% Yoel Shkolnisky, September 2013.

clear;

%% Generate simulated projections
% Generate 200 simulated projections of size 65x65.
% For simplicity, the projections are centered.
n=65;
K=200;
SNR=1/4;
%SNR=1000; % No noise
[projs,noisy_projs,~,rots_ref]=cryo_gen_projections(n,K,SNR);
viewstack(noisy_projs,5,5) % Show some noisy projections

masked_projs=mask_fuzzy(noisy_projs,33); % Applly circular mask
viewstack(masked_projs,5,5) % Show some noisy projections

% Compute polar Fourier transform, using radial resolution n_r and angular
% resolution n_theta. n_theta is the same as above.
n_theta=360;
n_r=65;
[npf,sampling_freqs]=cryo_pft(masked_projs,n_r,n_theta,'single');  % take Fourier transform of projections   

% Find common lines from projections
max_shift=0;
shift_step=1;
clstack = commonlines_gaussian(npf,max_shift,shift_step);

% Find reference common lines and compare
[ref_clstack,~]=clmatrix_cheat(rots_ref,n_theta);
prop=comparecl( clstack, ref_clstack, n_theta, 10 );
fprintf('Percentage of correct common lines: %f%%\n\n',prop*100);

%% Estimate orientations using sychronization.
 
S=cryo_syncmatrix_vote(clstack,n_theta);
[rotations,diff,mse]=cryo_syncrotations(S,rots_ref);
fprintf('MSE of the estimated rotations: %f\n\n',mse);
%fprintf('MSE of the estimated rotations: %f\n\n', ...
%   check_MSE(rotations,rots_ref));


%% 3D inversion
[ v, v_b, kernel ,err, iter, flag] = recon3d_firm( noisy_projs,...
rotations,[], 1e-6, 30, zeros(65,65,65));

assert(norm(imag(v(:)))/norm(v(:))<1.0e-3);
v=real(v);
WriteMRC(v,1,'example2.mrc'); % Output density map reconstructed from projections.
fprintf('Saved example2.mrc\n');
