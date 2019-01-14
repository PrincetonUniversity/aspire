
function [corrs_out]=clmatrix_ML_D2_par(pf,cls,doFilter,max_shift,shift_step)

ngpu=4;
n_proj=size(pf,3);
[idx_map]=lin2sub_map(n_proj);
corrs_out_par=cell(4,1);
npairs=nchoosek(n_proj,2);
% delete(gcp('nocreate'));
% gpuDevice([]);
% parpool('local',2);

if ~exist('max_shift','var')
    max_shift=0; % Maximal shift between common-lines in pixels. The 
                  % shift  is from -max_shift to max_shift. 
end

if ~exist('shift_step','var')
    shift_step=1.0; % Resolution of shift estimation in pixels.
end

tic
parfor i=1:4
    idx_map_loc=idx_map;    
%     if i==1
%         npairs1=floor(npairs/2);
%         [~,~,corrs_out_par{i}]=...
%             cryo_clmatrix_ML_gpu(pf,cls,n_proj,ngpu,idx_map_loc(1:npairs1,:),max_shift,shift_step,doFilter);
%     else
%         npairs1=floor(npairs/2);
%         [~,~,corrs_out_par{i}]=...
%             cryo_clmatrix_ML_gpu(pf,cls,n_proj,ngpu,idx_map_loc(npairs1+1:end,:),max_shift,shift_step,doFilter);
%     end
npairs1=floor(npairs/4);
if i<4
    [~,~,corrs_out_par{i}]=...
        cryo_clmatrix_ML_gpu(pf,cls,n_proj,ngpu,idx_map_loc((i-1)*npairs1+1:i*npairs1,:),...
        max_shift,shift_step,doFilter);
else
    [~,~,corrs_out_par{i}]=...
        cryo_clmatrix_ML_gpu(pf,cls,n_proj,ngpu,idx_map_loc((i-1)*npairs1+1:end,:),...
        max_shift,shift_step,doFilter);
end
end

corrs_out=corrs_out_par{1};
for i=2:4
    corrs_out.corrs=[corrs_out.corrs;corrs_out_par{i}.corrs];
    corrs_out.corrs_idx=[corrs_out.corrs_idx;corrs_out_par{i}.corrs_idx];
end
toc


