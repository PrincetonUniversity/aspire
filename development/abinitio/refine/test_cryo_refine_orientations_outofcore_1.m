Nprojs=100;
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
[Rs,shifts]=cryo_orient_projections_gpu(projshifted,voldata.volref,-1,trueRs,1,0);
t_orient=toc(t_orient);
fprintf('Assigning orientations took %5.1f seconds\n',t_orient);

% Refine orientations in-core 
initstate;
t_refined=tic;
[R_refined1,shifts_refined1,errs1]=cryo_refine_orientations(...
    projshifted,0,voldata.volref,Rs,shifts,1,-1,trueRs,true_shifts);
t_refined=toc(t_refined);
fprintf('Refining orientations %5.1f seconds\n',t_refined);

% Refine orientations out-of-core 
projs_fname=tempmrcsname;
imstackwriter=imagestackWriter(projs_fname,Nprojs);
imstackwriter.append(projshifted);
imstackwriter.close;

initstate;
t_refined=tic;
[R_refined2,shifts_refined2,errs2]=cryo_refine_orientations_outofcore(...
    projs_fname,0,voldata.volref,Rs,shifts,1,-1,trueRs,true_shifts);
t_refined=toc(t_refined);
fprintf('Refining orientations %5.1f seconds\n',t_refined);

% Print results
fprintf('Different should be 0. Results should match to the bit.\n');
fprintf('Difference in rotations = %e\n',norm(R_refined1(:)-R_refined2(:))/norm(R_refined1(:)));
fprintf('Difference in shifts = %e\n',norm(shifts_refined1(:)-shifts_refined2(:))/norm(shifts_refined1(:)));
fprintf('Difference in errors = %d\n',norm(errs1(:)-errs2(:))/norm(errs1(:)));

delete(projs_fname);
