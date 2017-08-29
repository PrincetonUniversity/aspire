% Test the function XXX
%
% Generate projections of a volume and estimate their orientations using
% the reference volume. Reconstruct the volume from the projections and
% their estimated orientations. Compare the reconstructed volume to the
% reference volume.
%
% Yoel Shkolnisky, June 2016.
clear;
initstate;
Nprojs=2000;
rots = rand_rots(Nprojs);  % Generate Nprojs projections to orient.

log_message('Loading volume');
voldata=load('cleanrib');

log_message('Generating %d clean projections',Nprojs);
projs=cryo_project(voldata.volref,rots);
projs=permute(projs,[2,1,3]);

log_message('Adding shifts');
[projshifted,ref_shifts]=cryo_addshifts(projs,[],2,1);

snr=1/8;
log_message('Adding noise. snr=%d',snr);
projshifted=cryo_addnoise(projshifted,snr,'gaussian');

% Invert rotations
trueRs=zeros(3,3,Nprojs);
for k=1:Nprojs
    trueRs(:,:,k)=rots(:,:,k).';
end

log_message('Orienting projecctions according to reference volume');
tic;
[Rest1,dx1]=cryo_orient_projections_gpu(projshifted,voldata.volref,[],trueRs,1,0);
t1_gpu=toc;

log_message('Results of the orientation procedure:')
log_message('Timing = %5.2f seconds',t1_gpu);

rot1_L2_error=norm(Rest1(:)-trueRs(:))/norm(trueRs(:));
log_message('L2 error in rotations estimation = %e',rot1_L2_error);

rot1_Killing_error=diag(dist_between_rot(Rest1,trueRs))/pi*180; %In degrees
log_message('Killing error in rotations estimation (in degrees)');
log_message('\t Max  = %5.3f',max(rot1_Killing_error));
log_message('\t mean = %5.3f',mean(rot1_Killing_error));
log_message('\t std  = %5.3f',std(rot1_Killing_error));
log_message('\t med  = %5.3f',median(rot1_Killing_error));

shifts1_L2_error=norm(dx1.'-ref_shifts)/norm(ref_shifts);
log_message('L2 error in shifts estimation= %e',shifts1_L2_error);
log_message('Max shift error in integral pixels (in each coordinate) = (%d,%d)',...
    max(round(ref_shifts)-round(dx1')));

log_message('Reconstructing from the projections and their etimated orientation parameters');
n=size(projshifted,1);
[ v1, ~, ~ ,~, ~, ~] = recon3d_firm( projshifted,Rest1,-dx1.', 1e-8, 100, zeros(n,n,n));
ii1=norm(imag(v1(:)))/norm(v1(:));
log_message('Relative norm of imaginary components = %e',ii1);
v1=real(v1);

log_message('Comparing reconstructed volume to reference volume');
plotFSC(voldata.volref,v1);

% Optimize rotations
rots_ref=rand_rots(Nprojs);  % Generate Nprojs projections to orient.
volref=voldata.volref;
projs_ref=cryo_project(volref,rots_ref);
projs_ref=permute(projs_ref,[2,1,3]);
Rrefs=zeros(3,3,Nprojs);
for k=1:Nprojs
    Rrefs(:,:,k)=rots_ref(:,:,k).';
end

L=360;
szvol=size(voldata.volref);
n_r=ceil(szvol(1)/2);
projs_ref_hat=cryo_pft(projs_ref,n_r,L,'single');
proj_hat=cryo_pft(projshifted,n_r,L,'single');

R_refined=zeros(size(Rest1));
estdx_refined=zeros(size(dx1));
errs=zeros(4,Nprojs);
parfor k=1:Nprojs
    estdx=dx1(:,k); 
    Rest=Rest1(:,:,k);
    [R_refined(:,:,k),estdx_refined(:,k),optout]=optimize_orientation(proj_hat(:,:,k),Rest,projs_ref_hat,Rrefs,L,estdx);
    
    rot_L2_error_before_refinement=norm(Rest1(:,:,k)-trueRs(:,:,k),'fro')/norm(trueRs(:,:,k),'fro');
    rot_L2_error_after_refinement=norm(R_refined(:,:,k)-trueRs(:,:,k),'fro')/norm(trueRs(:,:,k),'fro');

    shifts_L2_error_before_refinement=norm(dx1(:,k)-(ref_shifts(k,:)).');
    shifts_L2_error_after_refinement=norm(estdx_refined(:,k)-(ref_shifts(k,:)).');

    log_message('k=%d/%d',k,Nprojs);
    log_message('\t Rotation error: before=%e \t after=%e',rot_L2_error_before_refinement,rot_L2_error_after_refinement);
    log_message('\t Shift error:    before=%e \t after=%e',shifts_L2_error_before_refinement,shifts_L2_error_after_refinement);
    
    errs(:,k)=[shifts_L2_error_before_refinement;...
        shifts_L2_error_after_refinement;...
        shifts_L2_error_before_refinement;...
        shifts_L2_error_after_refinement];

end

[ v2, ~, ~ ,~, ~, ~] = recon3d_firm( projshifted,R_refined,-estdx_refined.', 1e-8, 100, zeros(n,n,n));
ii1=norm(imag(v2(:)))/norm(v2(:));
log_message('Relative norm of imaginary components = %e',ii1);
v2=real(v2);

log_message('Comparing reconstructed volume to reference volume');
plotFSC(voldata.volref,v2);

plotFSC2(voldata.volref,v1,voldata.volref,v2)