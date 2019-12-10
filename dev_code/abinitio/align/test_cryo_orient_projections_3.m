% Test the function cryo_orient_projections_gpu using a density map from
% EMDB.
%
% This script tests if we can match projections at one resolution to a
% density map at another resolution.  Currently this does
% not work. See DDD below.
%
% Yoel Shkolnisky, June 2016.

clear;
mapname=cryo_fetch_emdID(2275);
map=ReadMRC(mapname);
map=single(map);
%map=map.*fuzzymask(size(map),3,100,10);
%map=cryo_downsample(map,[65,65,65]);

% Generate projections at full resolution.
Nprojs=20;
initstate;
rots = rand_rots(Nprojs);  % Generate Nprojs projections to orient.
log_message('Generating %d projections to orient',Nprojs);
projs=cryo_project(map,rots);
projs=permute(projs,[2,1,3]);
log_message('Adding shifts');
[projshifted,ref_shifts]=cryo_addshifts(projs,[],2,1);
%[projshifted,ref_shifts]=cryo_addshifts(projs,[],0,1);

snr=1/8;
log_message('Adding noise wit snr=%d',snr);
projshifted=cryo_addnoise(projshifted,snr,'gaussian');

% Save trues rotations of the proejctions
% Invert rotations
trueRs=zeros(3,3,Nprojs);
for k=1:Nprojs
    trueRs(:,:,k)=rots(:,:,k).';
end


% Generate reference downsampled volume.
%map_downsampled=cryo_downsample(map,[129 129 129]);
% DDD
map_downsampled=map;


% Orient high-resolution projections
t_gpu=tic;
[Rest_gpu,dx_gpu]=cryo_orient_projections_gpu(projshifted,map_downsampled,[],trueRs,1,0);
%[Rest_gpu,dx_gpu]=cryo_orient_projections_gpu(projshifted,map,Nrefs,trueRs,0,0);
t_gpu=toc(t_gpu);

log_message('Timing = %5.2f seconds\n',t_gpu);
rot_L2_error=norm(Rest_gpu(:)-trueRs(:))/norm(trueRs(:));
log_message('L2 error in rotations estimation = %e\n',rot_L2_error);

rot_Killing_error=diag(dist_between_rot(Rest_gpu,trueRs))/pi*180; %In degrees
log_message('Killing error in rotations estimation (in degrees)\n');
log_message('\t Max  = %5.3f\n',max(rot_Killing_error));
log_message('\t mean = %5.3f\n',mean(rot_Killing_error));
log_message('\t std  = %5.3f\n',std(rot_Killing_error));
log_message('\t med  = %5.3f\n',median(rot_Killing_error));

shifts_L2_error=norm(dx_gpu.'-ref_shifts)/norm(ref_shifts);
log_message('L2 error in shifts estimation with preprocessing= %e\n',shifts_L2_error);
log_message('Max shift error in integral pixels (in each coordinate)= (%d,%d)\n',...
    max(round(ref_shifts)-round(dx_gpu.')));

