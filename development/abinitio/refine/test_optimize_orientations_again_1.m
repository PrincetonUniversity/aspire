Nprojs=10;
rots = rand_rots(Nprojs);  % Generate Nprojs projections to orient.
voldata=load('cleanrib');
projs=cryo_project(voldata.volref,rots);
projs=permute(projs,[2,1,3]);
[projshifted,true_shifts]=cryo_addshifts(projs,[],5,2);
true_shifts=true_shifts.';
snr=100000;
projshifted=cryo_addnoise(projshifted,snr,'gaussian');

% Invert rotations
trueRs=zeros(3,3,Nprojs);
for k=1:Nprojs
    trueRs(:,:,k)=rots(:,:,k).';
end

L = size(projs, 1);
[estR2,estdx2,optout]=optimize_orientations_again(projshifted,trueRs,true_shifts,L,trueRs,true_shifts);

rot_diff=norm(estR2(:)-trueRs(:))/norm(trueRs(:));
fprintf('Rotations difference from reference = %e\n',rot_diff);

shifts_diff=norm(estdx2(:)-true_shifts(:))/norm(true_shifts(:));
fprintf('Shifts difference from reference = %e\n',shifts_diff);
