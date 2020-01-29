
function [corrs_out]=clmatrix_ML_D2_scls_par(pf,cls,doFilter,max_shift,shift_step,scls_data,gpuIdx)
%% Default parameters.
if ~exist('max_shift','var')
    max_shift=0; % Maximal shift between common-lines in pixels. The 
                  % shift  is from -max_shift to max_shift. 
end

if ~exist('shift_step','var')
    shift_step=1.0; % Resolution of shift estimation in pixels.
end

%% Detect number of iterations needed for gpu's. 
%  Compute availible memory for use. 
scores=single(scls_data.scls_scores);
n_shifts=ceil(2*max_shift/shift_step+1);
nGpu = length(gpuIdx);
gpuMem = zeros(1,nGpu);
for j = 1:nGpu
   gpu_j = gpuDevice(gpuIdx(j));
   gpuMem(j) = gpu_j.AvailableMemory;
end
avlMem = min(gpuMem);
reqMem = 4*(3*length(cls)+3*length(cls)/4+numel(scores)+n_shifts*size(pf,2)^2);
nIter = ceil(reqMem/((0.95*avlMem))); % Each common line index on GPU is a single float, thus 4 
                                                  % bytes are required. We only fill pack 98%
                                                  % of memory on each GPU, and distribute
                                                  % operations on all avalible gpu's. 
nIter = max(nIter,nGpu);                                                  

%% Initialize parameters.
n_proj=size(pf,3);
npairs=nchoosek(n_proj,2);
corrs_out_par=cell(nIter,1);
tic

% Distribute ML computations between all GPU's by breaking the list of all
% common lines candidates (tuples) according to availible memory on arrays.
M=length(cls)/4;
n_per_gpu_iter=floor(M/nIter);
g_part=zeros(nIter,2);
g_part(1,:)=[1,n_per_gpu_iter];
for i=2:nIter-1
    g_part(i,:)=[g_part(i-1,2)+1,...
        g_part(i-1,2)+n_per_gpu_iter];
end
if nIter>1
    g_part(nIter,:)=[g_part(nIter-1,2)+1,M];
end
g_offset=g_part(:,1)-1;
l=g_part(:,2)-g_part(:,1)+1;

% Map the self common line scores of each 2 candidate rotations R_i,R_j to
% the respective realtive rotation candidate R_i^TR_J
nlookup1=length(scls_data.scls_lookup1)/3;
oct1_ij_map=scls_data.oct1_ij_map;
oct1_ij_map=[oct1_ij_map;oct1_ij_map(:,[2,1])];
oct2_ij_map=scls_data.oct2_ij_map;
oct2_ij_map(:,2)=oct2_ij_map(:,2)+nlookup1;
oct2_ij_map=[oct2_ij_map;oct2_ij_map(:,[2,1])];
ij_map=[oct1_ij_map;oct2_ij_map];
clear oct1_ij_map oct2_ij_map
scls_idx=zeros(l(1),2,nIter-1,'single');
for i=1:nIter-1
    scls_idx(:,:,i)=single(ij_map(g_part(i,1):g_part(i,2),:));
end
scls_idx_last=single(ij_map(g_part(nIter,1):g_part(nIter,2),:));

g_part(:,2)=g_part(:,2)*4;
g_part(:,1)=g_part(:,1)*4-3;
cl_idx=zeros(l(1)*4,nIter-1,'single');
for i=1:nIter-1
    cl_idx(:,i)=single(cls(g_part(i,1):g_part(i,2)));
end
cl_idx_last=single(cls(g_part(nIter,1):g_part(nIter,2)));

%% Run parallel code
log_message('Estimating D2 relative rotations');

delete(gcp('nocreate'));
parpool('local', nGpu);
spmd 
    gpuDevice(gpuIdx(labindex));
end

parfor i=1:nIter
    cl_idx_loc=cl_idx;
    scls_idx_loc=scls_idx;
    if i==nIter
        [corrs_out_par{i}]=...
            cryo_clmatrix_ML_gpu_scls(pf,cl_idx_last,max_shift,shift_step,...
            doFilter,scls_idx_last,scores);
    else
        [corrs_out_par{i}]=...
            cryo_clmatrix_ML_gpu_scls(pf,cl_idx_loc(:,i),max_shift,shift_step,...
            doFilter,scls_idx_loc(:,:,i),scores);
    end
end
toc
delete(gcp('nocreate'));

%Take max corr (score) over all nIter iterations
scores_agg=zeros(npairs,nIter,'single');
scores_idx=zeros(npairs,nIter);
for i=1:nIter
    scores_agg(:,i)=corrs_out_par{i}.corrs; 
    scores_idx(:,i)=corrs_out_par{i}.corrs_idx+g_offset(i);
end

[scores,best_iter]=max(scores_agg,[],2);
best_scores_offset=(1:npairs)'+(best_iter-1)*npairs;
tmp=scores_idx;
best_scores_idx=tmp(best_scores_offset);

corrs_out.corrs=scores;
corrs_out.corrs_idx=best_scores_idx;



