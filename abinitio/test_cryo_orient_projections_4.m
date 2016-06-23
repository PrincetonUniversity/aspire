% Test the function cryo_orient_projections_gpu.
%
% Generate projections of a volume and estimate their orientations using
% the reference volume. Reconstruct the volume from the projections and
% their estimated orientations. Compare the reconstructed volume to the
% reference volume.
%
% Yoel Shkolnisky, June 2016.

Nprojs=5000;
q=qrand(Nprojs);  % Generate Nprojs projections to orient.

log_message('Loading volume');
voldata=load('cleanrib');

log_message('Generating %d clean projections',Nprojs);
projs=cryo_project(voldata.volref,q);
projs=permute(projs,[2,1,3]);

log_message('Adding shifts');
[projshifted,ref_shifts]=cryo_addshifts(projs,[],2,1);

snr=1000;
log_message('Adding noise. snr=%d',snr);
projshifted=cryo_addnoise(projshifted,snr,'gaussian');

% Convert quaternions to rotations
trueRs=zeros(3,3,Nprojs);
for k=1:Nprojs
    trueRs(:,:,k)=(q_to_rot(q(:,k))).';
end

% Run with preprocessing
log_message('Orienting projecctions according to reference volume');
tic;
[Rest1_gpu,dx1_gpu]=cryo_orient_projections_gpu(projshifted,voldata.volref,[],trueRs,1,0);
t1_gpu=toc;

log_message('Results of the orientation procedure:')
fprintf('Timing = %5.2f seconds\n',t1_gpu);

rot1_L2_error=norm(Rest1_gpu(:)-trueRs(:))/norm(trueRs(:));
fprintf('L2 error in rotations estimation with preprocessing= %e\n',rot1_L2_error);

rot1_Killing_error=diag(dist_between_rot(Rest1_gpu,trueRs))/pi*180; %In degrees
fprintf('Killing error in rotations estimation (in degrees) with preprocessing\n');
fprintf('\t Max  = %5.3f\n',max(rot1_Killing_error));
fprintf('\t mean = %5.3f\n',mean(rot1_Killing_error));
fprintf('\t std  = %5.3f\n',std(rot1_Killing_error));
fprintf('\t med  = %5.3f\n',median(rot1_Killing_error));

shifts1_L2_error=norm(dx1_gpu.'+ref_shifts)/norm(ref_shifts);
fprintf('L2 error in shifts estimation with preprocessing= %e\n',shifts1_L2_error);
fprintf('Max shift error in integral pixels (in each coordinate) with preprocessing= (%d,%d)\n',...
    max(round(ref_shifts)+round(dx1_gpu')));

log_message('Reconstructing from the projections and their etimated orientation parameters');
n=size(projshifted,1);
[ v1, ~, ~ ,~, ~, ~] = recon3d_firm( projshifted,Rest1_gpu,dx1_gpu.', 1e-6, 30, zeros(n,n,n));
ii1=norm(imag(v1(:)))/norm(v1(:));
log_message('Relative norm of imaginary components = %e\n',ii1);
v1=real(v1);

log_message('Comparing reconstructed volume to reference volume');
plotFSC(voldata.volref,v1);
