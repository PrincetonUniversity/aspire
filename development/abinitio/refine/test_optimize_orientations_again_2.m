% Test the function optimize_orientations_again

Nprojs=100;
initstate;
rots = rand_rots(Nprojs);  % Generate Nprojs projections to orient.
voldata=load('cleanrib');
projs=cryo_project(voldata.volref,rots);
projs=permute(projs,[2,1,3]);
[projshifted,true_shifts]=cryo_addshifts(projs,[],2,1);
true_shifts=true_shifts.';
snr=1;
projshifted=cryo_addnoise(projshifted,snr,'gaussian');

% Invert rotations
trueRs=zeros(3,3,Nprojs);
for k=1:Nprojs
    trueRs(:,:,k)=rots(:,:,k).';
end

% Estimate rotations of the projections
t_orient=tic;
[Rs,shifts]=cryo_orient_projections(projshifted,voldata.volref,-1,trueRs,1,0);
t_orient=toc(t_orient);
log_message('Assigning orientations took %5.1f seconds',t_orient);


rot_diff=norm(Rs(:)-trueRs(:))/norm(trueRs(:));
log_message('Rotations difference from reference = %e',rot_diff);

shifts_diff=norm(shifts(:)-true_shifts(:))/norm(true_shifts(:));
log_message('Shifts difference from reference = %e',shifts_diff);

L=360;
[estR2,estdx2,optout]=optimize_orientations_again(projshifted,Rs,shifts,L,trueRs,true_shifts);

rot_diff=norm(estR2(:)-trueRs(:))/norm(trueRs(:));
log_message('Rotations difference from reference = %e',rot_diff);

shifts_diff=norm(estdx2(:)-true_shifts(:))/norm(true_shifts(:));
log_message('Shifts difference from reference = %e',shifts_diff);
