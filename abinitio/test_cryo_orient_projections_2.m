% Compare the performance of cryo_orient_projections_gpu with and without
% preprocessing of the projections.
%
% Note that preprocessing does not seem to help, so think about it some
% more.
%
% Yoel Shkolnisky, June 2016.


% Generate simulated data - proejctions and reference volume
Nprojs=10;
q=qrand(Nprojs);  % Generate Nprojs projections to orient.
voldata=load('cleanrib');
projs=cryo_project(voldata.volref,q);
projs=permute(projs,[2,1,3]);
[projshifted,ref_shifts]=cryo_addshifts(projs,[],2,1);
snr=1/20;
projshifted=cryo_addnoise(projshifted,snr,'gaussian');

% Convert quaternions to rotations
trueRs=zeros(3,3,Nprojs);
for k=1:Nprojs
    trueRs(:,:,k)=(q_to_rot(q(:,k))).';
end

Nrefs=100;
% Run with preprocessing
t1_gpu=tic;
[Rest1_gpu,dx1_gpu]=cryo_orient_projections_gpu(projshifted,voldata.volref,[],trueRs,1,1);
t1_gpu=toc(t1_gpu);

% Run without preprocessing
t2_gpu=tic;
[Rest2_gpu,dx2_gpu]=cryo_orient_projections_gpu(projshifted,voldata.volref,[],trueRs,1,0);
t2_gpu=toc(t2_gpu);


% Comprare results
fprintf('\n\n\n');
fprintf('Test results\n');
fprintf('************\n');

% rot_diff=norm(Rest1_gpu(:)-Rest2_gpu(:))/norm(Rest1_gpu(:));
% fprintf('Rotations difference = %e\n',rot_diff);
% 
% shifts_diff=norm(dx1_gpu(:)-dx2_gpu(:))/norm(dx1_gpu(:));
% fprintf('Shifts difference = %e\n',shifts_diff);

fprintf('Timing with preprocessing = %5.2f seconds\n',t1_gpu);
fprintf('Timing without preprocessing = %5.2f seconds\n',t2_gpu);

rot1_L2_error=norm(Rest1_gpu(:)-trueRs(:))/norm(trueRs(:));
fprintf('L2 error in rotations estimation with preprocessing= %e\n',rot1_L2_error);
rot2_L2_error=norm(Rest2_gpu(:)-trueRs(:))/norm(trueRs(:));
fprintf('L2 error in rotations estimation without preprocessing= %e\n',rot2_L2_error);

rot1_Killing_error=diag(dist_between_rot(Rest1_gpu,trueRs))/pi*180; %In degrees
fprintf('Killing error in rotations estimation (in degrees) with preprocessing\n');
fprintf('\t Max  = %5.3f\n',max(rot1_Killing_error));
fprintf('\t mean = %5.3f\n',mean(rot1_Killing_error));
fprintf('\t std  = %5.3f\n',std(rot1_Killing_error));
fprintf('\t med  = %5.3f\n',median(rot1_Killing_error));

rot2_Killing_error=diag(dist_between_rot(Rest2_gpu,trueRs))/pi*180; %In degrees
fprintf('Killing error in rotations estimation (in degrees) without preprocessing\n');
fprintf('\t Max  = %5.3f\n',max(rot2_Killing_error));
fprintf('\t mean = %5.3f\n',mean(rot2_Killing_error));
fprintf('\t std  = %5.3f\n',std(rot2_Killing_error));
fprintf('\t med  = %5.3f\n',median(rot2_Killing_error));

shifts1_L2_error=norm(dx1_gpu.'-ref_shifts)/norm(ref_shifts);
fprintf('L2 error in shifts estimation with preprocessing= %e\n',shifts1_L2_error);
fprintf('Max shift error in integral pixels (in each coordinate) with preprocessing= (%d,%d)\n',...
    max(round(ref_shifts)-round(dx1_gpu')));

shifts2_L2_error=norm(dx2_gpu.'-ref_shifts)/norm(ref_shifts);
fprintf('L2 error in shifts estimation without preprocessing= %e\n',shifts2_L2_error);
fprintf('Max shift error in integral pixels (in each coordinate) without preprocessing= (%d,%d)\n',...
    max(round(ref_shifts)-round(dx2_gpu.')));
