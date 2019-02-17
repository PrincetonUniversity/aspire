
function [corrs_out]=clmatrix_ML_Dn_scls_par(pf,cls,doFilter,max_shift,shift_step,scls_data,n)

n_proj=size(pf,3);
npairs=nchoosek(n_proj,2);
ngpu=16;
corrs_out_par=cell(ngpu,1);
%whichScore=3;

if ~exist('max_shift','var')
    max_shift=0; % Maximal shift between common-lines in pixels. The 
                  % shift  is from -max_shift to max_shift. 
end

if ~exist('shift_step','var')
    shift_step=1.0; % Resolution of shift estimation in pixels.
end
tic
M=length(cls)/(2*n);
n_per_gpu_iter=floor(M/ngpu);
g_part=zeros(ngpu,2);
g_part(1,:)=[1,n_per_gpu_iter];
for i=2:ngpu-1
    g_part(i,:)=[g_part(i-1,2)+1,...
        g_part(i-1,2)+n_per_gpu_iter];
end
if ngpu>1
    g_part(ngpu,:)=[g_part(ngpu-1,2)+1,M];
end
g_offset=g_part(:,1)-1;
l=g_part(:,2)-g_part(:,1)+1;

%Handle self common line scores
nlookup1=length(scls_data.scls_lookup1)/(2*n);
scores=single(scls_data.scls_scores);
oct1_ij_map=scls_data.oct1_ij_map;
oct1_ij_map=[oct1_ij_map;oct1_ij_map(:,[2,1])];
oct2_ij_map=scls_data.oct2_ij_map;
oct2_ij_map(:,2)=oct2_ij_map(:,2)+nlookup1;
oct2_ij_map=[oct2_ij_map;oct2_ij_map(:,[2,1])];
ij_map=[oct1_ij_map;oct2_ij_map];
clear oct1_ij_map oct2_ij_map
scls_idx=zeros(l(1),2,ngpu-1,'single');
for i=1:ngpu-1
    scls_idx(:,:,i)=single(ij_map(g_part(i,1):g_part(i,2),:));
end
scls_idx_last=single(ij_map(g_part(ngpu,1):g_part(ngpu,2),:));

%g_part(:,2)=g_part(:,2)*4;
%g_part(:,1)=g_part(:,1)*4-3;
g_part(:,2)=g_part(:,2)*(2*n);
g_part(:,1)=g_part(:,1)*(2*n)-(2*n-1);
%g_part=[0,0;g_part];
cl_idx=zeros(l(1)*(2*n),ngpu-1,'single');
for i=1:ngpu-1
    cl_idx(:,i)=single(cls(g_part(i,1):g_part(i,2)));
end
cl_idx_last=single(cls(g_part(ngpu,1):g_part(ngpu,2)));

parfor i=1:ngpu
    cl_idx_loc=cl_idx;
    scls_idx_loc=scls_idx;
    if i==ngpu
        [~,~,corrs_out_par{i}]=...
            cryo_clmatrix_ML_gpu_scls_Dn(pf,cl_idx_last,max_shift,shift_step,...
            doFilter,scls_idx_last,scores,n);
    else
        [~,~,corrs_out_par{i}]=...
            cryo_clmatrix_ML_gpu_scls_Dn(pf,cl_idx_loc(:,i),max_shift,shift_step,...
            doFilter,scls_idx_loc(:,:,i),scores,n);
    end
end
toc
% corrs_out=[];
% return;

%Take max corr (score) over all ngpu iterations
scores_agg=zeros(npairs,ngpu,3,'single');
scores_idx=zeros(npairs,ngpu,3);
for i=1:ngpu
    scores_agg(:,i,:)=corrs_out_par{i}.corrs; 
    scores_idx(:,i,:)=corrs_out_par{i}.corrs_idx+g_offset(i);
end

scores=zeros(npairs,3,'single');
best_scores_idx=zeros(npairs,3);

for i=1:3
    [scores(:,i),best_iter]=max(scores_agg(:,:,i),[],2);
    best_scores_offset=(1:npairs)'+(best_iter-1)*npairs;
    tmp=scores_idx(:,:,i);
    best_scores_idx(:,i)=tmp(best_scores_offset);    
end
corrs_out.corrs=scores;
corrs_out.corrs_idx=best_scores_idx;

% corrs_out=corrs_out_par{1};
% corrs_out.corrs=[corrs_out.corrs;corrs_out_par{2}.corrs];
% corrs_out.corrs_idx=[corrs_out.corrs_idx;corrs_out_par{2}.corrs_idx];
% toc


