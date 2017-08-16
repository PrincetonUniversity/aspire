% Verify that the output of the functions
%       cryo_orient_projections_ref
%       cryo_orient_projections
%       cryo_orient_projections_gpu
% are the same.
% This is a sanity function to verify changes to any of these functions.
%
% Yoel Shkolnisky, June 2016.


% Generate simulated data - proejctions and reference volume
Nprojs=20;
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

Nrefs=10;
% Run reference code
initstate; % All three functions must generate exactly the same temporary  data.
t_ref=tic;
[Rest_ref,dx_ref]=cryo_orient_projections_ref(projshifted,voldata.volref,Nrefs,trueRs,1);
t_ref=toc(t_ref);

% Run CPU code
initstate;
t_cpu=tic;
[Rest_cpu,dx_cpu]=cryo_orient_projections(projshifted,voldata.volref,Nrefs,trueRs,1);
t_cpu=toc(t_cpu);

% Run GPU code
initstate;
t_gpu=tic;
[Rest_gpu,dx_gpu]=cryo_orient_projections_gpu(projshifted,voldata.volref,Nrefs,trueRs,1);
t_gpu=toc(t_gpu);

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
fprintf('GPU: Errors of the order of 10^-8 are allowed since tables are stored as single\n');


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


shifts_L2_error=norm(dx_gpu.'-ref_shifts)/norm(ref_shifts);
fprintf('L2 error in shifts estimation = %e\n',shifts_L2_error);
fprintf('Max shift error in integral pixels (in each coordinate) = (%d,%d)\n',...
    max(round(ref_shifts)-round(dx_gpu')));

% Check precomputed tables
% Compare refernce and CPU tables
tables_ref=load(fullfile(tempmrcdir,'cryo_orient_projections_tables_ref.mat'));
tables_cpu=load(fullfile(tempmrcdir,'cryo_orient_projections_tables_cpu.mat'));
assert(all(tables_ref.Cjk(:)==tables_cpu.Cjk(:)));
assert(all(tables_ref.Ckj(:)==tables_cpu.Ckj(:)));
assert(all(tables_ref.Mkj(:)==tables_cpu.Mkj(:)));
assert(all(tables_ref.qrefs(:)==tables_cpu.qrefs(:)));

% Compare refernce and GPU tables
tables_ref=load(fullfile(tempmrcdir,'cryo_orient_projections_tables_ref.mat'));
tables_gpu=load(fullfile(tempmrcdir,'cryo_orient_projections_tables_gpu.mat'));
assert(all(tables_ref.Cjk(:)==tables_gpu.Cjk(:)));
assert(all(tables_ref.Ckj(:)==tables_gpu.Ckj(:)));
assert(all(tables_ref.Mkj(:)==tables_gpu.Mkj(:)));
assert(all(tables_ref.qrefs(:)==tables_gpu.qrefs(:)));
log_message('Test OK');
