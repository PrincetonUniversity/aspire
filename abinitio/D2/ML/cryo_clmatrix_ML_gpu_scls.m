function [corrs_data]=...
    cryo_clmatrix_ML_gpu_scls(pf,cl_matrix,max_shift,shift_step,doFilter,scls_matrix,scores)
%
%
%   Run common lines Maximum likelihood procedure for a D2 molecule, to
%   find the set of rotations Ri^TgkRj, k=1,2,3,4 for each pair of images i
%   and j. 
%
% Input parameters:
%   pf       3D array where each image pf(:,:,k) corresponds to the Fourier
%            transform of projection k.
%   max_shift       Maximal 1D shift (in pixels)  to search between
%       common-lines. Default: 15.
%   shift_step      Resolution of shift estimation in pixels. Note that
%        shift_step can be any positive real number. Default: 1.
%   cl_matrix       a vector with linear 
% Returned variables:
%
%   corrs_data      For each pair of images i and j return the best
%   estimate (with highest common line correlations score) for the relative 
%   rotations Ri^Tg_kRj (k=1,2,3,4). 


%PRECISION='single';
initstate;
T=size(pf,2);
if mod(T,2)~=0
    error('n_theta must be even');
end

% pf is of size n_rxn_theta. Convert pf into an array of size
% (2xn_r-1)xn_theta, that is, take then entire ray through the origin, but
% thake the angles only up PI.
% This seems redundant: The original projections are real, and thus
% each ray is conjugate symmetric. We therefore gain nothing by taking
% longer correlations (of length 2*n_r-1 instead of n_r), as the two halfs
% are exactly the same. Taking shorter correlation would speed the
% computation by a factor of two.
pf=[flipdim(pf(2:end,T/2+1:end,:),1) ; pf(:,1:T/2,:) ];
n_theta=size(pf,2);
n_proj=size(pf,3);

%% Check input parameters and set debug flags.
if ~exist('max_shift','var');
    max_shift=15; % Maximal shift between common-lines in pixels. The 
                  % shift  is from -max_shift to max_shift. 
end

if ~exist('shift_step','var');
    shift_step=1.0; % Resolution of shift estimation in pixels.
end
n_shifts=ceil(2*max_shift/shift_step+1); % Number of shifts to try.

%% Allocate output variables
clstack=zeros(n_proj,n_proj);      % Common lines-matrix.
corrstack=zeros(n_proj,n_proj);    % Correlation coefficient for each common-line.                               
%% Search for common lines between pairs of projections

% Construct filter to apply to each Fourier ray.                   
rmax=(size(pf,1)-1)/2;    
rk=-rmax:rmax; rk=rk(:);
H=sqrt(abs(rk)).*exp(-rk.^2/(2*(rmax/4).^2)); 
H=repmat(H(:),1,n_theta);  % Filter for common-line detection.

% Bandpass filter and normalize each ray of each projection.
% XXX We do not override pf since it is used to debugging plots below. Once
% XXX these debugging plots are removed, replace pf3 by pf. This will save
% XXX a lot of memory. 
pf3=pf;
if doFilter
    for k=1:n_proj
        proj=pf(:,:,k);
        proj=proj.*H;
        proj(rmax:rmax+2,:)=0;
        proj=cryo_raynormalize(proj);
        pf3(:,:,k)=proj;
    end
else
    for k=1:n_proj
        proj=pf(:,:,k);
        proj(rmax:rmax+2,:)=0;
        proj=cryo_raynormalize(proj);
        pf3(:,:,k)=proj;
    end
end

rk2=rk(1:rmax);
% Prepare the shift_phases
shift_phases=zeros(rmax,n_shifts);
for shiftidx=1:n_shifts
    shift=-max_shift+(shiftidx-1)*shift_step;
    shift_phases(:,shiftidx)=exp(-2*pi*sqrt(-1).*rk2.*shift./(2*rmax+1));
end

%% Run ML in parallel
npairs=nchoosek(n_proj,2);
corrs_idx_out=zeros(npairs,1);
corrs_out=zeros(npairs,1,'single');
ij_idx=0;
l=gpuArray(single(length(cl_matrix)/4));
cl_matrix_g=gpuArray(cl_matrix);
scls_matrix_i=gpuArray(scls_matrix(:,1));
scls_matrix_j=gpuArray(scls_matrix(:,2));
scores_g=gpuArray(scores);

for k1=1:n_proj-1

    proj1=pf3(:,:,k1);
    P1=proj1(1:rmax,:);  % Take half ray plus the DC
    P1_flipped=conj(P1);
    
    % Instead of processing the shifts in a loop, one shift at a time, stack
    % all shifts of a single projection into one array and compute all
    % correlations for all shifts at once.
    % The variable P1_stack stacks all shifted versions of the image k1.
    P1_stack=zeros(size(P1,1),2*size(P1,2)*n_shifts);
    for k=1:n_shifts
        P1_stack(:,(2*k-2)*n_theta+1:(2*k-1)*n_theta)=bsxfun(@times,P1,shift_phases(:,k));      
        P1_stack(:,(2*k-1)*n_theta+1:(2*k)*n_theta)=bsxfun(@times,P1_flipped,shift_phases(:,k));
    end

    g_P1_stack=gpuArray(single(P1_stack));

    % Make sure the DC component is zero. This is assumed  below in
    % computing correlations.
    if norm(proj1(rmax+1,:))>1.0e-13
        error('DC component of projection is not zero');
    end
    
    scores_i=scores_g(:,k1);
    
    for k2=k1+1:n_proj
        ij_idx=ij_idx+1;                    
        proj2=pf3(:,:,k2); % proj1 and proj2 are both normalized to unit norm.
        P2=proj2(1:rmax,:);
        g_P2=gpuArray(single(P2));  
        
        if norm(proj2(rmax+1,:))>1.0e-13
            error('DC component of projection is not zero');
        end               
        
        g_C=2*real(g_P1_stack'*g_P2);
        g_C=reshape(g_C,2*n_theta,n_shifts,n_theta);
        g_C=permute(g_C,[1,3,2]);
        g_C=squeeze(max(g_C,[],3));
        %g_C=0.5*(g_C+1); %map to probabilities
        
        tmp=g_C(cl_matrix_g);
        tmp=reshape(tmp,4,l);
        prod_corrs=prod(tmp,1);
        clearvars tmp
        
        %Incorporate scores of individual rotations calculated from self-cls
        scores_j=scores_g(:,k2);
        scores_ij=scores_i(scls_matrix_i);
        scores_ij=scores_ij.*scores_j(scls_matrix_j);
        
        prod_corrs=prod_corrs'.*scores_ij;
        [max_corr,max_idx]=max(prod_corrs);    
        corrs_idx_out(ij_idx)=gather(max_idx);
        corrs_out(ij_idx)=gather(max_corr);        
    end
end        
corrs_data=struct('corrs_idx',corrs_idx_out,'corrs',corrs_out);


