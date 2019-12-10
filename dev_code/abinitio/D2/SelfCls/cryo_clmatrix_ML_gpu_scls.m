function [clstack,corrstack,corrs_data]=...
    cryo_clmatrix_ML_gpu_scls(pf,cl_matrix,max_shift,shift_step,doFilter,scls_matrix,scores)
%
%
%   Generate common-lines matrix for the Fourier stack pf.
%   This function is identical to cryo_clmatrix but takes advantage of GPU.
%   See cryo_clmatrix for details.
%   This version (using single precision GPU arithmetic) is 2.4 times
%   faster than cryo_clmatrix (at version 990 at SVN).
%   Note that when using single precision there are slight discrepencies
%   between this function and cryo_clmatrix. When changing all GPU code to
%   double precision, the result of this function is identical to
%   cryo_clmatrix.
%
% Input parameters:
%   pf       3D array where each image pf(:,:,k) corresponds to the Fourier
%            transform of projection k.
%   NK       For each projection find its common-lines with NK other
%            projections. If NK is less than the total number a projection,
%            a random subset of NK projections is used. Default: n_proj. 
%   max_shift       Maximal 1D shift (in pixels)  to search between
%       common-lines. Default: 15.
%   shift_step      Resolution of shift estimation in pixels. Note that
%        shift_step can be any positive real number. Default: 1.
%
% Returned variables:
%
%   clstack     Common lines matrix. (k1,k2) and (k2,k1) contain the index
%       of the common line of projections k1 and k2. (k1,k2)  contains the
%       index of the common line in the projection k1. (k2,k1) contains the
%       index of the common line in k2. 
%   corrstack   The correlation of the common line between projections k1
%       and k2. Since corrstack is symmetric, it contain entries only above
%       the diagonal. corrstack(k1,k2) measures how ''common'' is the between
%       projections k1 and k2. Small value means high-similariry.
%   shift_equations  System of equations for the 2D shifts of the
%       projections. This is a sparse system with 2*n_proj+1 columns. The
%       first 2*n_proj columns correspond to the unknown shifts dx,dy of
%       each projection. The last column is the right-hand-side of the
%       system, which is the relative shift of each pair of common-lines.
%   shift_equations_map   2D array of size n_proj by n_proj. Entry (k1,k2)
%       is the index of the equation (row number) in the array
%       "shift_equations" that corresponds to the common line between
%       projections k1 and k2. shift_map is non-zero only for k1<k2. 


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
corrs_idx_out=zeros(npairs,3);
corrs_out=zeros(npairs,3,'single');
ij_idx=0;
l=gpuArray(single(length(cl_matrix)/4));
cl_matrix_g=gpuArray(cl_matrix);
scls_matrix_i=gpuArray(scls_matrix(:,1));
scls_matrix_j=gpuArray(scls_matrix(:,2));
scores_g=gpuArray(scores);
%deb_arr=zeros(npairs,3);

for k1=1:n_proj-1
    %k1=ij_map(p_idx,1);
    %waitbar(ij_idx/npairs);

    proj1=pf3(:,:,k1);
    P1=proj1(1:rmax,:);  % Take half ray plus the DC
    P1_flipped=conj(P1);
    
    % Instead of processing the shifts in a loop, one shift at a time, stack
    % all shifts of a single projection into one array and compute all
    % correlations for all shifts at once.
    % The variable P1_stack stacks all shifted versions of the image k1.
    P1_stack=zeros(size(P1,1),2*size(P1,2)*n_shifts);
    %P1_flipped_stack=zeros(size(P1,1),size(P1,2)*n_shifts);
    for k=1:n_shifts
        P1_stack(:,(2*k-2)*n_theta+1:(2*k-1)*n_theta)=bsxfun(@times,P1,shift_phases(:,k));
        %P1_flipped_stack(:,(k-1)*n_theta+1:k*n_theta)=bsxfun(@times,P1_flipped,shift_phases(:,k));
        P1_stack(:,(2*k-1)*n_theta+1:(2*k)*n_theta)=bsxfun(@times,P1_flipped,shift_phases(:,k));
    end

    g_P1_stack=gpuArray(single(P1_stack));
%     if strcmpi(PRECISION,'single')
%         g_P1_stack=gpuArray(single(P1_stack));        
%         %g_P1_flipped_stack=gpuArray(single(P1_flipped_stack));
%     else
%         g_P1_stack=gpuArray(P1_stack);
%         %g_P1_flipped_stack=gpuArray(P1_flipped_stack);
%     end

    % Make sure the DC component is zero. This is assumed  below in
    % computing correlations.
    if norm(proj1(rmax+1,:))>1.0e-13
        error('DC component of projection is not zero');
    end
    
    scores_i=scores_g(:,k1);
    
    for k2=k1+1:n_proj
    %k2=ij_map(p_idx,2);
        ij_idx=ij_idx+1;                    
        proj2=pf3(:,:,k2); % proj1 and proj2 are both normalized to unit norm.
        P2=proj2(1:rmax,:);
        g_P2=gpuArray(single(P2));  
        
%         if strcmpi(PRECISION,'single')
%             g_P2=gpuArray(single(P2));            
%         else
%             g_P2=gpuArray(P2);
%         end
%         
        if norm(proj2(rmax+1,:))>1.0e-13
            error('DC component of projection is not zero');
        end               
        
%Optimized GPU version:
        g_C=2*real(g_P1_stack'*g_P2);
        g_C=reshape(g_C,2*n_theta,n_shifts,n_theta);
        g_C=permute(g_C,[1,3,2]);
        g_C=squeeze(max(g_C,[],3));
        g_C=0.5*(g_C+1); %map to probabilities
        
        tmp=g_C(cl_matrix_g);
        tmp=reshape(tmp,4,l);
        %tmp=log10(tmp);
        %sum_corrs=sum(tmp,1);
        %sum_corrs=sum_corrs'/4;
        prod_corrs=prod(tmp,1);
        %prod_corrs=sum(tmp,1);
        %prod_corrs=geomean(tmp+1,1)-1;
        clearvars tmp
        
        %Add scores of individual rotations calculated from self-cls
        scores_j=scores_g(:,k2);
        scores_ij=scores_i(scls_matrix_i);
        scores_ij=scores_ij.*scores_j(scls_matrix_j);
        %scores_ij=log10(scores_ij)+log10(scores_j(scls_matrix_j));
        
%         [max_corr,max_idx]=max(sum_corrs);  
%         corrs_idx_out(ij_idx,1)=gather(max_idx);
%         corrs_out(ij_idx,1)=gather(max_corr);      
        
        %[max_corr,max_idx]=max(scores_ij);  
        %corrs_idx_out(ij_idx,2)=gather(max_idx);
        %corrs_out(ij_idx,2)=gather(max_corr); 
%         sum_corrs=sum_corrs+scores_i(scls_matrix_i);
%         sum_corrs=sum_corrs+scores_j(scls_matrix_j);    

        %sum_corrs=sum_corrs+scores_ij;
        %sum_corrs=sum_corrs.*scores_ij;
        %[max_corr,max_idx]=max(sum_corrs);    
        %corrs_out(ij_idx,4)=gather(scores_ij(max_idx));
        %corrs_idx_out(ij_idx,3)=gather(max_idx);
        %corrs_out(ij_idx,3)=gather(max_corr);  
        %clear sum_corrs
        
        prod_corrs=prod_corrs'.*scores_ij;
        %prod_corrs=prod_corrs'+scores_ij;
        [max_corr,max_idx]=max(prod_corrs);    
        %prod_corrs2=gather(prod_corrs);
        %[max_corr2,max_idx2]=max(prod_corrs2);
        %max_idx3=gather(max_idx);
        corrs_idx_out(ij_idx,2)=gather(max_idx);
        corrs_out(ij_idx,2)=gather(max_corr);
       
        
        %DEBUG
        %deb_arr(ij_idx,2)=gather(ceil(scls_matrix_j(max_idx)/72));
        %deb_arr(ij_idx,1)=gather(ceil(scls_matrix_i(max_idx)/72));
        %deb_arr(ij_idx,3)=double(gather(max_corr));
        
        sval=gather(max_corr);
        sidx=gather(max_idx);
        
        clstack(k1,k2)=sidx;
        corrstack(k1,k2)=sval;
%         corrs_idx_out(ij_idx)=sidx;
%         corrs_out(ij_idx)=sval;          
    end
end
%corrs_out=corrs_out(:,3);
%corrs_idx_out=corrs_idx_out(:,3);
corrstack(corrstack~=0)=1-corrstack(corrstack~=0);            
corrs_data=struct('corrs_idx',corrs_idx_out,'corrs',corrs_out);


