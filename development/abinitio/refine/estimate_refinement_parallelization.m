function [numworkers,t1,t2]=estimate_refinement_parallelization(n)
%
% ESTIMATE_REFINEMENT_PARALLELIZATION   Estimate parallelization for max
%           GPU utilization
%
% [numworkers,t1,t2]=estimate_refinement_parallelization(n)
%   Estimate the number of concurrent processes that should be used by
%   cryo_orient_projections_gpu to maximize GPU utilization.
%   Returns the optimal number of workers, the timing t1 using 1 worker (no
%   palallelization) and t2 which is an array with two columns where the
%   first column is the number of workers used and the second is the timing
%   for that number of workers. The retuned numworkers is the one for which
%   the speedup relative to t1 is maximal.
%   
% Yoel Shkolnisky, July 2016.

Nprojs=200;
q=qrand(Nprojs);  % Generate Nprojs projections to orient.

log_message('Loading volume');

% Maps (EMD codes) and their size
%   2660    360
%   6474    320
%   2275    240
%   3418    160

if n>360
    error('n larger than 360 is currently not supported.');
end

mapname=cryo_fetch_emdID(2660);
vol=ReadMRC(mapname);
vol=cryo_downsample(vol,[n n n],0);


log_message('Generating %d clean projections of size %dx%d',Nprojs,size(vol,1),size(vol,2));
projs=cryo_project(vol,q);
projs=permute(projs,[2,1,3]);

log_message('Adding shifts');
[projshifted,~]=cryo_addshifts(projs,[],2,1);

snr=500;
log_message('Adding noise. snr=%d',snr);
projshifted=cryo_addnoise(projshifted,snr,'gaussian');

% Convert quaternions to rotations
trueRs=zeros(3,3,Nprojs);
for k=1:Nprojs
    trueRs(:,:,k)=(q_to_rot(q(:,k))).';
end

log_message('Generating reference timing using single worker.');
t1=tic;
[Rest1_gpu,dx1_gpu]=cryo_orient_projections_gpu(projshifted,vol,[],trueRs,1,0);
t1=toc(t1);

candidate_numworkers=[2 4 6 8 12];
%candidate_numworkers=[12];
t2=zeros(numel(candidate_numworkers),2);
for k=1:numel(candidate_numworkers)
    log_message('Timing using %d workers.',candidate_numworkers(k));
    t2(k,1)=candidate_numworkers(k);
    t2(k,2)=-1;
    try
        poolreopen(candidate_numworkers(k));
        t2_gpu=tic;
        [Rest2_gpu,dx2_gpu]=cryo_orient_projections_gpu_2(projshifted,vol,[],trueRs,1,0,candidate_numworkers(k));
        t2_gpu=toc(t2_gpu);
        assert(norm(Rest1_gpu(:)-Rest2_gpu(:))==0);
        assert(norm(dx1_gpu(:)-dx2_gpu(:))==0);        
        t2(k,2)=t2_gpu;
    catch E
        log_message('Error occured for numworkers=%d:\n%s',candidate_numworkers(k),E.identifier);
    end
end

t=t1./t2(:,2);
numworkers=candidate_numworkers(t==max(t) & t>0);