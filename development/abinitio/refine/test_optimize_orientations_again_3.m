Nprojs=100;
q=qrand(Nprojs);  % Generate Nprojs projections to orient.
voldata=load('cleanrib');
projs=cryo_project(voldata.volref,q);
projs=permute(projs,[2,1,3]);
[projshifted,true_shifts]=cryo_addshifts(projs,[],2,1);
true_shifts=true_shifts.';
snr=1/8;
projshifted=cryo_addnoise(projshifted,snr,'gaussian');

% Convert quaternions to rotations
trueRs=zeros(3,3,Nprojs);
for k=1:Nprojs
    trueRs(:,:,k)=(q_to_rot(q(:,k))).';
end
L=360;
n_r=ceil(size(projshifted,1)/2);
pf=cryo_pft(projshifted,n_r,L,'single');
clmatrix=cryo_clmatrix(pf,Nprojs,1,6,1);
S=cryo_syncmatrix_vote(clmatrix,L);
Rest=cryo_syncrotations(S);
dxest=cryo_estimate_shifts(pf,Rest,6,1);

figure(1)
dir1=Q2S2(q,L);
dir2=R2S2(Rest,L);
check_orientations(dir1,dir2);

[regrot,mse]=register_rotations(Rest,trueRs);
                                
rot_diff=norm(regrot(:)-trueRs(:))/norm(trueRs(:));
fprintf('Rotations difference from reference = %e\n',rot_diff);

L=360;
[Rest2,dxest2,optout]=optimize_orientations_again(projshifted,Rest,dxest,L,trueRs,true_shifts);

figure(2)
dir1=Q2S2(q,L);
dir2=R2S2(Rest2,L);
check_orientations(dir1,dir2);


rot_diff=norm(Rest2(:)-trueRs(:))/norm(trueRs(:));
fprintf('Rotations difference from reference = %e\n',rot_diff);

shifts_diff=norm(dxest2(:)-true_shifts(:))/norm(true_shifts(:));
fprintf('Shifts difference from reference = %e\n',shifts_diff);
