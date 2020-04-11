function [corrs_out]=...
    cryo_clmatrix_scl_ML(pf,max_shift,shift_step,doFilter,scls_data)
%
% Input parameters:
%   pf       3D array where each image pf(:,:,k) corresponds to the Fourier
%            transform of projection k.
%   max_shift       Maximal 1D shift (in pixels)  to search between
%       common-lines. Default: 15.
%   shift_step      Resolution of shift estimation in pixels. Note that
%        shift_step can be any positive real number. Default: 1.
%
% Output parameters: 
%   


%PRECISION='single';
5;
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
if ~exist('max_shift','var')
    max_shift=15; % Maximal shift between common-lines in pixels. The 
                  % shift  is from -max_shift to max_shift. 
end

if ~exist('shift_step','var');
    shift_step=1.0; % Resolution of shift estimation in pixels.
end
n_shifts=ceil(2*max_shift/shift_step+1); % Number of shifts to try.
     
%% Search for common lines between pairs of projections

% Construct filter to apply to each Fourier ray.                   
rmax=(size(pf,1)-1)/2;    
rk=-rmax:rmax; rk=rk(:);
H=sqrt(abs(rk)).*exp(-rk.^2/(2*(rmax/4).^2)); 
H=repmat(H(:),1,n_theta);  % Filter for common-line detection.

if doFilter
    for k=1:n_proj
        proj=pf(:,:,k);
        proj=proj.*H;
        proj(rmax:rmax+2,:)=0;
        proj=cryo_raynormalize(proj);
        pf(:,:,k)=proj;
    end
else
    for k=1:n_proj
        proj=pf(:,:,k);
        proj(rmax:rmax+2,:)=0;
        proj=cryo_raynormalize(proj);
        pf(:,:,k)=proj;
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
cl_matrix=[scls_data.scls_lookup1;scls_data.scls_lookup2];
M=length(cl_matrix)/3;
corrs_out=zeros(M,n_proj,'single');
cl_idx=reshape(cl_matrix,3,M);

%Get linear indices to use with corrs matrix
nonEq_idx=scls_data.nonEq_idx;
nonTv_eq_idx=scls_data.nonTv_eq_idx;
nonEq_lin_idx=cl_idx(:,nonEq_idx);
nonEq_lin_idx=nonEq_lin_idx(:);
nNonEq=length(nonEq_idx);
%tv_idx=scls_data.tv_idx;

%Get lists of candidate indices for equators
eq_idx_lists=scls_data.eq_idx_lists;
ntheta=scls_data.ntheta;
n_eq=scls_data.n_eq;
% dtheta=round(180*scls_data.dtheta/pi);
% n_tv=scls_data.n_tv;
% tv_zero_angles=scls_data.tv_zero_angles;

for k1=1:n_proj
    
    proj1=pf(:,:,k1);
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
    
    proj2=pf(:,:,k1); % proj1 and proj2 are both normalized to unit norm.
    P2=proj2(1:rmax,:);
    g_P2=gpuArray(single(P2));
    
    if norm(proj2(rmax+1,:))>1.0e-13
        error('DC component of projection is not zero');
    end
    
    % Compute correlations between all common lines
    g_C=2*real(g_P1_stack'*g_P2);
    g_C=reshape(g_C,2*n_theta,n_shifts,n_theta);
    g_C=permute(g_C,[1,3,2]);
    g_C=squeeze(max(g_C,[],3));
    
    %map correaltions to probabilities (in the spirit of Maximum Likelihood)
    g_C=0.5*(g_C+1);
    
    %Compute equator measures
    [eq_measures,normals_corrs_max]=allEqMeasures(g_C,360);
    eq_measures=eq_measures.*normals_corrs_max; %This can be done inside allEqmeasure. 
    %Handle cases: 
    %1. Non equators - just take product of probabilities
    prod_corrs=prod(reshape(g_C(nonEq_lin_idx),3,nNonEq),1);
    corrs_out(nonEq_idx,k1)=gather(prod_corrs);
    %2. Non tv equators - adjust scores by eq measures
    for eq_idx=1:n_eq
        for i=1:ntheta           
           % Take the correlations for the self common line candidate of the
           % "equator rotation" #eq_idx with respect to image k1, and 
           % multiply by all scores from the function eq_measures (see 
           % documentation inside the function ). Then take maximum over
           % all the scores. 
           true_scls_corrs=g_C(eq_idx_lists{i,eq_idx,1});                       
           scls_cand_idx=eq_idx_lists{i,eq_idx,2};
           eq_measures_i=eq_measures(scls_cand_idx);            
           measures_agg=true_scls_corrs'*eq_measures_i;
           j=nonTv_eq_idx(eq_idx);
           corrs_out((j-1)*ntheta+i,k1)=gather(max(measures_agg(:))); %DO: Check if this is correct indexing
           %This line can be made more efficient by preallocating and
           %avoiding gather calls
        end      
    end  
    
    
    %3. Top view (tv) equators DO: Need to compute only one tv_ind, since
    %they depend only on in plane rotation all tv are the same    
%     tv_measures=gather(tv_measures);
%     for tv_ind=1:n_tv
%         j=tv_idx(tv_ind);
%         %corrs_out((j-1)*ntheta+(1:ntheta),k1)=tv_measures;
%         idx1=mod((tv_zero_angles(tv_ind,1))+360-(0:ntheta-1)*dtheta,360);
%         idx1(idx1==0)=360;
%         idx2=mod((tv_zero_angles(tv_ind,2))+360-(0:ntheta-1)*dtheta,360);
%         idx2(idx2==0)=360;        
%         corrs_out((j-1)*ntheta+(1:ntheta),k1)=...
%             max([tv_measures(idx1);tv_measures(idx2)],[],1);
%         for i=1:ntheta
%             theta=(i-1)*dtheta;
%             j=tv_idx(tv_ind);
%             %corrs_out((j-1)*ntheta+i,k1)=gather(measureTv(g_C,theta,L));
%             corrs_out((j-1)*ntheta+i,k1)=gather(tv_measures(theta));
%         end
%    end
end


