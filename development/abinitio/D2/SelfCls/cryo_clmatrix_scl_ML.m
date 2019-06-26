function [corrs_out]=...
    cryo_clmatrix_scl_ML(pf,cl_matrix,max_shift,shift_step,doFilter,scls_data)
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
M=length(cl_matrix)/3;
corrs_out=zeros(M,n_proj,'single');
cl_idx=reshape(cl_matrix,3,M);
%cl_idx=gpuArray(single(reshape(cl_matrix,3,M)));

%Get linear indices to use with corrs matrix
nonEq_idx=scls_data.nonEq_idx;
nonTv_eq_idx=scls_data.nonTv_eq_idx;
tv_idx=scls_data.tv_idx;
nonEq_lin_idx=cl_idx(:,nonEq_idx);
nonEq_lin_idx=nonEq_lin_idx(:);
nNonEq=length(nonEq_idx);

%Get lists of candidate indices for equators
eq_idx_lists=scls_data.eq_idx_lists;
ntheta=scls_data.ntheta;
dtheta=round(180*scls_data.dtheta/pi);
% L=scls_data.L;
n_eq=scls_data.n_eq;
n_tv=scls_data.n_tv;
tv_zero_angles=scls_data.tv_zero_angles;

for k1=1:n_proj
    
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
    
    % Make sure the DC component is zero. This is assumed  below in
    % computing correlations.
    if norm(proj1(rmax+1,:))>1.0e-13
        error('DC component of projection is not zero');
    end
    
    proj2=pf3(:,:,k1); % proj1 and proj2 are both normalized to unit norm.
    P2=proj2(1:rmax,:);
    g_P2=gpuArray(single(P2));
    
    if norm(proj2(rmax+1,:))>1.0e-13
        error('DC component of projection is not zero');
    end
    
    g_C=2*real(g_P1_stack'*g_P2);
    g_C=reshape(g_C,2*n_theta,n_shifts,n_theta);
    g_C=permute(g_C,[1,3,2]);
    g_C=squeeze(max(g_C,[],3));
    
    %map to probabilities
    g_C=0.5*(g_C+1);
    
    %Compute equator measures
    [eq_measures,normals_corrs_max,scls_corrs,tv_measures]=allEqMeasures(g_C,360);
    eq_measures=eq_measures.*normals_corrs_max;
    %Handle cases: 
    %1. Non equators
    prod_corrs=prod(reshape(g_C(nonEq_lin_idx),3,nNonEq),1);
    corrs_out(nonEq_idx,k1)=gather(prod_corrs);
    %2. Non tv equators    
    for eq_idx=1:n_eq
        for i=1:ntheta           
            true_scls_corrs=g_C(eq_idx_lists{i,eq_idx,1});                       
            scls_cand_idx=eq_idx_lists{i,eq_idx,2};
            eq_measures_i=eq_measures(scls_cand_idx);
%            l=length(scls_cand_idx);
%            eqm=zeros(1,l);
%             for k=1:l
%                 m=measureEq2(g_C,scls_cand_idx(l),L);
%                 eqm(k)=m(1)*m(3);
%             end
%            tmp=true_scls_corrs'*eqm;  
           j=nonTv_eq_idx(eq_idx);
           measures_agg=true_scls_corrs'*eq_measures_i;
           corrs_out((j-1)*ntheta+i,k1)=gather(max(measures_agg(:))); %DO: Check if this is correct indexing
            %This line can be made more efficient by preallocating and
            %avoiding gather calls
        end      
    end  
    
    
    %3. Top view (tv) equators DO: Need to compute only one tv_ind, since
    %they depend only on in plane rotation all tv are the same    
    tv_measures=gather(tv_measures);
    for tv_ind=1:n_tv
        j=tv_idx(tv_ind);
        %corrs_out((j-1)*ntheta+(1:ntheta),k1)=tv_measures;
        idx1=mod((tv_zero_angles(tv_ind,1))+360-(0:ntheta-1)*dtheta,360);
        idx1(idx1==0)=360;
        idx2=mod((tv_zero_angles(tv_ind,2))+360-(0:ntheta-1)*dtheta,360);
        idx2(idx2==0)=360;        
        corrs_out((j-1)*ntheta+(1:ntheta),k1)=...
            max([tv_measures(idx1);tv_measures(idx2)],[],1);
%         for i=1:ntheta
%             theta=(i-1)*dtheta;
%             j=tv_idx(tv_ind);
%             %corrs_out((j-1)*ntheta+i,k1)=gather(measureTv(g_C,theta,L));
%             corrs_out((j-1)*ntheta+i,k1)=gather(tv_measures(theta));
%         end
    end
        
    %mean_corrs=g_C(cl_idx);
    %mean_corrs=reshape(mean_corrs,3,M);
    %mean_corrs=sum(mean_corrs,1)/3;
    %corrs_out(:,k1)=gather(mean_corrs);
    %[m,I]=max(mean_corrs);
end


