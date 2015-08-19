% Generate simulated data - proejctions and reference volume
Nprojs=10;
q=qrand(Nprojs);  % Generate Nprojs projections to orient.
voldata=load('cleanrib');
projs=cryo_project(voldata.volref,q);
projs=permute(projs,[2,1,3]);
[projshifted,ref_shifts]=cryo_addshifts(projs,[],2,1);
snr=1000;
projshifted=cryo_addnoise(projshifted,snr,'gaussian');

% Convert quaternions to rotations
trueRs=zeros(3,3,Nprojs);
for k=1:Nprojs
    trueRs(:,:,k)=(q_to_rot(q(:,k))).';
end

% Run reference code
tic;
[Rest_ref,dx_ref]=cryo_orient_projections_ref(projshifted,voldata.volref,10,trueRs,0);
t_ref=toc;

% Run CPU code
tic;
[Rest_cpu,dx_cpu]=cryo_orient_projections(projshifted,voldata.volref,10,trueRs,0);
t_cpu=toc;

% Run GPU code
tic;
[Rest_gpu,dx_gpu]=cryo_orient_projections_gpu(projshifted,voldata.volref,10,trueRs,0);
t_gpu=toc;

% Comprare results
fprintf('\n\n\n');
fprintf('Test results\n');
fprintf('************\n');

rot_diff=norm(Rest_cpu(:)-Rest_ref(:))/norm(Rest_ref(:));
fprintf('CPU: Rotations difference from reference = %e\n',rot_diff);

shifts_diff=norm(dx_cpu(:)-dx_ref(:))/norm(dx_ref(:));
fprintf('CPU: Shifts difference from reference = %e\n',shifts_diff);

rot_diff=norm(Rest_gpu(:)-Rest_ref(:))/norm(Rest_ref(:));
fprintf('GPU: Rotations difference from reference = %e\n',rot_diff);

shifts_diff=norm(dx_gpu(:)-dx_ref(:))/norm(dx_ref(:));
fprintf('GPU: Shifts difference from reference = %e\n',shifts_diff);


fprintf('Timing of reference code = %5.2f seconds\n',t_ref);
fprintf('Timing of CPU  code = %5.2f seconds\n',t_cpu);
fprintf('Timing of GPU  code = %5.2f seconds\n',t_gpu);

rot_L2_error=norm(Rest_gpu(:)-trueRs(:))/norm(trueRs(:));
fprintf('L2 error in rotations estimation = %e\n',rot_L2_error);

rot_Killing_error=diag(dist_between_rot(Rest_gpu,trueRs))/pi*180; %In degrees
fprintf('Killing error in rotations estimation (in degrees)\n');
fprintf('\t Max  = %5.3f\n',max(rot_Killing_error));
fprintf('\t mean = %5.3f\n',mean(rot_Killing_error));
fprintf('\t std  = %5.3f\n',std(rot_Killing_error));


shifts_L2_error=norm(dx_gpu.'+ref_shifts)/norm(ref_shifts);
fprintf('L2 error in shifts estimation = %e\n',shifts_L2_error);
fprintf('Max shift error in integral pixels (in each coordinate) = (%d,%d)\n',...
    max(round(ref_shifts)+round(dx_gpu')));
