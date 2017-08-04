
clear;
open_log(0);

N=100;
SNR=1;
n_theta=360;
max_shift=0;
shift_step=1;

[~,stack,~,rots_ref]=cryo_gen_projections(89,N,SNR,max_shift,shift_step);
[ref_clstack,~]=clmatrix_cheat_q(rot_to_q(rots_ref),n_theta);

% Normalize background noise to mean 0 and std 1
stack=cryo_normalize_background(stack);
n_r=ceil(size(stack,1)*0.5);

% Mask projections
K=size(stack,3);
mask_radius=round(size(stack,1)*0.45);
[masked_stack,~]=mask_fuzzy(stack,mask_radius);
figure(1);viewstack(masked_stack,5,5);

% Use correlation
[pf_corr,~]=cryo_pft(masked_stack,n_r,n_theta,'single');  % take Fourier transform of projections
clstack_corr = cryo_clmatrix_gpu(pf_corr,K,1,round(2*sqrt(2)*max_shift),shift_step);
prop_corr=comparecl( clstack_corr, ref_clstack, n_theta, 5 );
log_message('Correlation: %f%%',prop_corr*100);

% Estimate noise power spectrum
spnoise_pft= cryo_noise_estimation_pfs(stack,n_r,n_theta);
spnoise_pft(spnoise_pft<max(spnoise_pft/10))=1;
noise_cov=diag(spnoise_pft);

% Polar Fourier transform
[pf_ml,sampling_freqs]=cryo_pft(masked_stack,n_r,n_theta,'single');  % take Fourier transform of projections

% Normalize the polar Fourier transform.
for k=1:size(pf_ml,3);
    X=masked_stack(:,:,k);
    Y=pf_ml(:,:,k);
    cX=norm(X(:));
    cY=norm(Y(:));
    Y=(Y./cY).*cX;
    pf_ml(:,:,k)=Y;
end

open_log(0);
max_shift=0;
shift_step=1;
max_iterations=6;
M=26000;
KNN=200;
[clstack_ml,~,shift_equations,shift_equations_map,~]=...
    cryo_clmatrix_ml_gpu(pf_ml,M,KNN,noise_cov,max_iterations,...
    1,max_shift,shift_step);

prop_ml=comparecl( clstack_ml, ref_clstack, n_theta, 5 );
log_message('ML: %f%%',prop_ml*100);
