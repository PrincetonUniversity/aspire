% Compute the polar Fourier transform outside cryo_refine_orientations_outofcore
Nprojs=1000;
rots = rand_rots(Nprojs);  % Generate Nprojs projections to orient.
voldata=load('cleanrib');
projs=cryo_project(voldata.volref,rots);
projs=permute(projs,[2,1,3]);
[projshifted,true_shifts]=cryo_addshifts(projs,[],2,1);
true_shifts=true_shifts.';
snr=1000;
projshifted=cryo_addnoise(projshifted,snr,'gaussian');
projshifted=single(projshifted);

% Invert rotations
trueRs=zeros(3,3,Nprojs);
for k=1:Nprojs
    trueRs(:,:,k)=rots(:,:,k).';
end

% Estimate rotations of the projections
t_orient=tic;
[Rs,shifts]=cryo_orient_projections_gpu(projshifted,voldata.volref,-1,trueRs,0,0);
t_orient=toc(t_orient);
fprintf('Assigning orientations took %5.1f seconds\n',t_orient);

% Refine orientations out-of-core
projs_fname=tempmrcsname;
imstackwriter=imagestackWriter(projs_fname,Nprojs);
imstackwriter.append(projshifted);
imstackwriter.close;

initstate;
t_refined=tic;
% XXX Note that L=360 is set by the following function internally. It is
% better to pass it as a parameter to match the L value below.
[R_refined1,shifts_refined1,errs1]=cryo_refine_orientations_outofcore(...
    projs_fname,0,voldata.volref,Rs,shifts,1,-1,trueRs,true_shifts);
t_refined=toc(t_refined);
fprintf('Refining orientations %5.1f seconds\n',t_refined);


% Compute polar Fourier transform of the projecitons.
L=360;
n_r=ceil(size(voldata.volref,1)/2);
projs_hat_fname=tempmrcsname;
cryo_pft_outofcore(projs_fname,projs_hat_fname,n_r,L);
projs_hat_normalized_fname=tempmrcsname;
cryo_raynormalize_outofcore(projs_hat_fname,projs_hat_normalized_fname); 
    % Fourier transformed projections are assumed to be normalized.

initstate;
t_refined=tic;
[R_refined2,shifts_refined2,errs2]=cryo_refine_orientations_outofcore(...
    projs_hat_normalized_fname,1,voldata.volref,Rs,shifts,1,-1,trueRs,true_shifts);
t_refined=toc(t_refined);
fprintf('Refining orientations %5.1f seconds\n',t_refined);

% Print results
fprintf('Different should be 0. Results should match to the bit.\n');
fprintf('Difference in rotations = %e\n',norm(R_refined1(:)-R_refined2(:))/norm(R_refined1(:)));
fprintf('Difference in shifts = %e\n',norm(shifts_refined1(:)-shifts_refined2(:))/norm(shifts_refined1(:)));
fprintf('Difference in errors = %d\n',norm(errs1(:)-errs2(:))/norm(errs1(:)));

delete(projs_fname);
delete(projs_hat_fname);
delete(projs_hat_normalized_fname);
