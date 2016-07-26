
 INPUTDIR='/home/yoel/data/';%projs500.mrc
 OUTPUTDIR='/tmp/clml';

 open_log(0);

    %% abinitio reconstruction
N=500;
%fname='projs500.mrc';
fname='projs500_simulated.mrc';
average=ReadMRC(fullfile(INPUTDIR,fname));
stack=average(:,:,1:N); 
% Normalize background noise to mean 0 and std 1
stack=cryo_normalize_background(stack);
n_theta=360;
n_r=ceil(size(stack,1)*0.5);

clear average;

% Mask projections
K=size(stack,3);
mask_radius=round(size(stack,1)*0.45);

[masked_stack,~]=mask_fuzzy(stack,mask_radius);

figure(1);viewstack(masked_stack,5,5);

% Estimate noise power spectrum
spnoise_pft= cryo_noise_estimation_pfs(stack,n_r,n_theta);
spnoise_pft(spnoise_pft<max(spnoise_pft/10))=1;
noise_cov=diag(spnoise_pft);

% Polar Fourier transform
[pf,sampling_freqs]=cryo_pft(stack,n_r,n_theta,'single');  % take Fourier transform of projections

% Normalize the polar Fourier transform.
for k=1:size(pf,3);
    X=stack(:,:,k);
    Y=pf(:,:,k);
    cX=norm(X(:));
    cY=norm(Y(:));
    Y=(Y./cY).*cX;
    pf(:,:,k)=Y;
end

% tmppf=reshape(pf,[size(pf,1) size(pf,2)*size(pf,3)]);
% pfnorms=sqrt(sum(abs(tmppf).^2,1));
% pf=pf./max(pfnorms);

% pf=cryo_raynormalize(reshape(pf,n_r,n_theta*N));
% pf=reshape(pf,n_r,n_theta,N);
% Find common lines from projections using maximum
% likelihood
% XXX REPLACE HERE XXX
open_log(0);
max_shift=10;
shift_step=1;
max_iterations=6;
M=26000;
KNN=50;
[clstack,~,shift_equations,shift_equations_map,~]=...
    cryo_clmatrix_ml_gpu(pf,M,KNN,noise_cov,max_iterations,...
    1,max_shift,shift_step);
%[clstack,~,shift_equations,shift_equations_map]=...
%    cryo_clmatrix_gpu(pf,K,1,max_shift,shift_step);
%%
reloadname=sprintf('reload_real_data_eee');
save(fullfile(OUTPUTDIR,reloadname));

S=cryo_syncmatrix_vote(clstack,n_theta);
rotations=cryo_syncrotations(S);

clerr=syncconsistency(rotations,clstack,n_theta);
figure;hist(clerr(:,3),360);

save(fullfile(OUTPUTDIR,reloadname));

[est_shifts,~]=cryo_estimate_shifts(pf,rotations,max_shift,shift_step,10000,[],0);
save(fullfile(OUTPUTDIR,reloadname));

% Reconstruct
n=size(stack,1);
[ v1, v_b1, kernel1 ,err1, iter1, flag1] = recon3d_firm( stack,...
    rotations,-est_shifts, 1e-6, 30, zeros(n,n,n));
ii1=norm(imag(v1(:)))/norm(v1(:));
fprintf('Relative norm of imaginary components = %e\n',ii1);
v1=real(v1);
volname=sprintf('projs500.mrc');
WriteMRC(v1,1,fullfile(OUTPUTDIR,volname));