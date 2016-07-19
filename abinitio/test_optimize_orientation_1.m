clear;
Nprojs=10;
q_ref=qrand(Nprojs);  % Generate Nprojs projections to orient.
voldata=load('cleanrib');
volref=voldata.volref;
projs_ref=cryo_project(volref,q_ref);
projs_ref=permute(projs_ref,[2,1,3]);

q=qrand(1);
proj=cryo_project(volref,q);
proj=permute(proj,[2,1,3]);
[proj,ref_shifts]=cryo_addshifts(proj,[],2,1);
% snr=1000;
% projshifted=cryo_addnoise(projshifted,snr,'gaussian');

trueRs=zeros(3,3,Nprojs);
for k=1:Nprojs
    trueRs(:,:,k)=(q_to_rot(q_ref(:,k))).';
end

R=(q_to_rot(q)).';

L=360;
szvol=size(volref);
n_r=ceil(szvol(1)/2);
projs_ref_hat=cryo_pft(projs_ref,n_r,L,'single');
proj_hat=cryo_pft(proj,n_r,L,'single');
estdx=ref_shifts.';



[Rest,estdx,optout]=optimize_orientation(proj_hat,R,projs_ref_hat,trueRs,L,estdx);

rot_L2_error=norm(Rest(:)-R(:))/norm(R(:));
fprintf('L2 error in rotations estimation = %e\n',rot_L2_error);

fprintf('ref_shifts=(%d,%d)\n',ref_shifts);
fprintf('estdx=(%5.3e,%5.3e)\n',estdx);
shifts_L2_error=norm(estdx.'-ref_shifts)/norm(ref_shifts);
fprintf('L2 error in shifts estimation = %e\n',shifts_L2_error);
fprintf('Max shift error in integral pixels (in each coordinate) = (%d,%d)\n',...
    (round(ref_shifts)-round(estdx.')));

