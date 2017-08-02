Nprojs=10;
rots = rand_rots(Nprojs);  % Generate Nprojs projections to orient.
voldata=load('cleanrib');
projs=cryo_project(voldata.volref,rots);
projs=permute(projs,[2,1,3]);
[projshifted,true_shifts]=cryo_addshifts(projs,[],2,1);
true_shifts=true_shifts.';
snr=1000;
projshifted=cryo_addnoise(projshifted,snr,'gaussian');

% Invert rotations
trueRs=zeros(3,3,Nprojs);
for k=1:Nprojs
    trueRs(:,:,k)=rots(:,:,k)';
end

% Estimate rotations of the projections
t_orient=tic;
[Rs,shifts]=cryo_orient_projections(projshifted,voldata.volref,-1,trueRs,1,0);
t_orient=toc(t_orient);
fprintf('Assigning orientations took %5.1f seconds\n',t_orient);

% Refine orientations
t_refined=tic;
[R_refined,shifts_refined,errs]=cryo_refine_orientations(...
    projshifted,0,voldata.volref,Rs,shifts,1,-1,trueRs,true_shifts);
t_refined=toc(t_refined);
fprintf('Refining orientations %5.1f seconds\n',t_refined);

% Somehow show that the errors get smaller after refinement, e.g.,
median(errs)
