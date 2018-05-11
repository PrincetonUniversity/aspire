% Compare the output of cryo_orient_projections_gpu against
% cryo_orient_projections_gpu_outofcore.
%
% The output of both functions should be excatly the same (to the bit).
%
% Yoel Shkolnisky, April 2017.


%% Generate simulated data - proejctions and reference volume
Nprojs=1000;
rots = rand_rots(Nprojs);  % Generate Nprojs projections to orient.
voldata=load('cleanrib');
projs=cryo_project(voldata.volref,rots);
projs=permute(projs,[2,1,3]);
[projshifted,ref_shifts]=cryo_addshifts(projs,[],2,1);
snr=1000;
projshifted=cryo_addnoise(projshifted,snr,'gaussian');
projshifted=single(projshifted);

% Invert rotations
trueRs=zeros(3,3,Nprojs);
for k=1:Nprojs
    trueRs(:,:,k)=rots(:,:,k).';
end

Nrefs=10;

%% Run cryo_orient_projections_gpu

t_gpu=tic;
initstate;
[Rest_gpu,dx_gpu]=cryo_orient_projections_gpu(projshifted,voldata.volref,Nrefs,trueRs,1);
t_gpu=toc(t_gpu);

%% Run cryo_orient_projections_gpu_outofcore

% Write projections to MRC file
projs_fname='temp.mrcs';
projswriter=imagestackWriter(projs_fname,Nprojs);
projswriter.append(projshifted);
projswriter.close;

% Orient projections
t_gpu_outofcore=tic;
initstate;
[Rest_gpu_outofcore,dx_gpu_outofcore]=...
    cryo_orient_projections_gpu_outofcore(projs_fname,voldata.volref,Nrefs,trueRs,1);
t_gpu_outofcore=toc(t_gpu_outofcore);


% Comprare results
fprintf('\n\n\n');
fprintf('Test results\n');
fprintf('************\n');

rot_diff=norm(Rest_gpu(:)-Rest_gpu_outofcore(:))/norm(Rest_gpu(:));
fprintf('Rotations difference = %e\n',rot_diff);

shifts_diff=norm(dx_gpu(:)-dx_gpu_outofcore(:))/norm(dx_gpu(:));
fprintf('Shifts difference = %e\n',shifts_diff);
fprintf('GPU: Errors of the order of 10^-8 are allowed since tables are stored as single\n');

fprintf('Timing of GPU code = %5.2f seconds\n',t_gpu);
fprintf('Timing of OFC code = %5.2f seconds\n',t_gpu_outofcore);
