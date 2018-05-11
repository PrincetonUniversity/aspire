Nprojs=10;
initstate;
rots = rand_rots(Nprojs);  % Generate Nprojs projections to orient.
voldata=load('cleanrib');
projs=cryo_project(voldata.volref,rots);
projs=permute(projs,[2,1,3]);
[projshifted,ref_shifts]=cryo_addshifts(projs,[],2,1);
snr=1000;
projshifted=cryo_addnoise(projshifted,snr,'gaussian');

% Invert rotations
trueRs=zeros(3,3,Nprojs);
for k=1:Nprojs
    trueRs(:,:,k)=rots(:,:,k).';
end

Nrefs=10;
[Rest_gpu,dx_gpu]=cryo_orient_projections_gpu(projshifted,voldata.volref,Nrefs,trueRs,1);

rot_L2_error=norm(Rest_gpu(:)-trueRs(:))/norm(trueRs(:));
log_message('L2 error in rotations estimation = %e',rot_L2_error);
shifts_L2_error=norm(dx_gpu.'-ref_shifts)/norm(ref_shifts);
log_message('L2 error in shifts estimation = %e',shifts_L2_error);
log_message('Max shift error in integral pixels (in each coordinate) = (%d,%d)',...
    max(round(ref_shifts)-round(dx_gpu')));

rots_ref = rand_rots(Nprojs);  % Generate Nprojs projections to orient.
volref=voldata.volref;
projs_ref=cryo_project(volref,rots_ref);
projs_ref=permute(projs_ref,[2,1,3]);
Rrefs=zeros(3,3,Nprojs);
for k=1:Nprojs
    Rrefs(:,:,k)=rots_ref(:,:,k).';
end


L=360;
szvol=size(volref);
n_r=ceil(szvol(1)/2);
projs_ref_hat=cryo_pft(projs_ref,n_r,L,'single');
proj_hat=cryo_pft(projshifted,n_r,L,'single');

for k=1:Nprojs
    estdx=dx_gpu(:,k); 
    Rest=Rest_gpu(:,:,k);
    [R_refined,estdx_refined,optout]=optimize_orientation(proj_hat(:,:,k),Rest,projs_ref_hat,Rrefs,L,estdx);
    
    rot_L2_error_before_refinement=norm(Rest_gpu(:,:,k)-trueRs(:,:,k),'fro')/norm(trueRs(:,:,k),'fro');
    rot_L2_error_after_refinement=norm(R_refined-trueRs(:,:,k),'fro')/norm(trueRs(:,:,k),'fro');

    shifts_L2_error_before_refinement=norm(dx_gpu(:,k)-(ref_shifts(k,:)).');
    shifts_L2_error_after_refinement=norm(estdx_refined-(ref_shifts(k,:)).');

    log_message('k=%d/%d',k,Nprojs);
    log_message('\t Rotation error: before=%e \t after=%e',rot_L2_error_before_refinement,rot_L2_error_after_refinement);
    log_message('\t Shift error:    before=%e \t after=%e',shifts_L2_error_before_refinement,shifts_L2_error_after_refinement);

end
