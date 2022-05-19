function [clstack,corrstack]=...
    cryo_clmatrix_gpu_kernel_2(pf,NK,verbose,max_shift,shift_step,map_filter_radius)
%
%
%   Generate common-lines matrix for the Fourier stack pf.
%
% Input parameters:
%   pf       3D array where each image pf(:,:,k) corresponds to the Fourier
%            transform of projection k.
%   NK       For each projection find its common-lines with NK other
%            projections. If NK is less than the total number a projection,
%            a random subset of NK projections is used. Default: n_proj. 
%   verbose  Bitmask of debugging level.Bits:
%           0   silent
%           1   One line progress message (not written to log) (Default)
%   max_shift       Maximal 1D shift (in pixels)  to search between
%       common-lines. Default: 15.
%   shift_step      Resolution of shift estimation in pixels. Note that
%        shift_step can be any positive real number. Default: 1.
%   map_filter_radius      If nonzero, the common line between a pair
%       images is detected not by the pair of lines with the highest 
%       correlation, but rather the pair of lines that both them and their
%       sorroundings (map_filter_radius Fourier rays to each directions)
%       give the best match. The radius for comparison is
%       determined by the value of map_filter_radius (Default 0).
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
%
% Yoel Shkolnisky, May 2022

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
NK=n_proj;
if (nargin<2) || (NK==-1)
    NK=n_proj; % Number of common-line pairs to compute for each projection
end

if ~exist('verbose','var')
    verbose=1;
end

if ~exist('max_shift','var')
    max_shift=15; % Maximal shift between common-lines in pixels. The 
                  % shift  is from -max_shift to max_shift. 
end

if ~exist('shift_step','var')
    shift_step=1.0; % Resolution of shift estimation in pixels.
end
n_shifts=ceil(2*max_shift/shift_step+1); % Number of shifts to try.

if ~exist('map_filter_radius','var')
    map_filter_radius=0;
end


% Set flag for progress and debug messages
verbose_progress=0;

if bitand(verbose,1)
    verbose_progress=1;
end

if verbose~=0
    log_message('Verbose mode=%d',verbose);
end

%%

clstack=zeros(n_proj,n_proj);      % Common lines-matrix.
corrstack=zeros(n_proj,n_proj);    % Correlation coefficient for each common-line.

if verbose>0
    log_message('map_filter_radius = %d',map_filter_radius);
end

           
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
for k=1:n_proj
    proj=pf(:,:,k);
    proj=proj.*H;
    proj(rmax:rmax+2,:)=0;
    proj=cryo_raynormalize(proj);
    pf3(:,:,k)=proj;
end

kernel = parallel.gpu.CUDAKernel('kernel7.ptx','kernel7.cu');
kernel.ThreadBlockSize = [1, 1, 1];
kernel.GridSize = [1, 1];
g_clstack=gpuArray(zeros(n_proj,'single'));
g_corrstack=gpuArray(zeros(n_proj,'single'));

t_cpu=0;
t_gpu=0;

rk2=rk(1:rmax);
for k1=1:n_proj
   
    proj1=pf3(:,:,k1);
    P1=proj1(1:rmax,:);  % Take half ray plus the DC
    P1_flipped=conj(P1);
    
    % Make sure the DC component is zero. This is assumed  below in
    % computing correlations.
    if norm(proj1(rmax+1,:))>1.0e-13
        error('DC component of projection is not zero');
    end
    
    for k2=k1+1:n_proj 
        
        proj2=pf3(:,:,k2); % proj1 and proj2 are both normalized to unit norm.
        P2=proj2(1:rmax,:);
        
        if norm(proj2(rmax+1,:))>1.0e-13
            error('DC component of projection is not zero');
        end
            
         t1=clock;
         g_P1=gpuArray(single(P1));
         g_P2=gpuArray(single(P2));
 
         [g_clstack,g_corrstack]=feval(kernel,k1,k2,n_proj,g_P1,g_P2,rmax,n_theta,...
            max_shift, shift_step,g_clstack,g_corrstack);
         t2=clock;
         t_gpu=t_gpu+etime(t2,t1);
    
    %kernel = parallel.gpu.CUDAKernel('kernel6.ptx','kernel6.cu');
    %g_shift_phases=gpuArray(complex(zeros(rmax,n_shifts,'single')));

    %kernel.ThreadBlockSize = [1, 1, 1];
    %kernel.GridSize = [1, 1];

    %g_shift_phases=feval(kernel,max_shift, shift_step, rmax, g_shift_phases);
    
        t1=clock;
        % Find the shift that gives best correlation.
        for shiftidx=1:n_shifts
            shift=-max_shift+(shiftidx-1)*shift_step;
            shift_phases=exp(-2*pi*sqrt(-1).*rk2.*shift./(2*rmax+1)); 
            %disp(norm(g_shift_phases(:,shiftidx)-shift_phases)/norm(shift_phases));
            %shift_phases=repmat(shift_phases,1,n_theta);
            
            % No need to renormalize proj1_shifted and
            % proj1_shifted_flipped since multiplication by phases
            % does not change the norm, and proj1 is already normalized.
            %P1_shifted=P1.*shift_phases;
            %P1_shifted_flipped=P1_flipped.*shift_phases;
                                    
            P1_shifted=bsxfun(@times,P1,shift_phases);
            P1_shifted_flipped=bsxfun(@times,P1_flipped,shift_phases);
            
            % Compute correlations in the positive r direction           
            C1=2*real(P1_shifted'*P2);
            
            % Compute correlations in the negative r direction
            C2=2*real(P1_shifted_flipped'*P2);
                                   
            C = [C1,C2];
   
            if map_filter_radius > 0
                C = cryo_average_clmap(C, map_filter_radius);
            end
            
            [sval,sidx]=max(C(:));            
           
            if sval>corrstack(k1,k2)
                [cl1,cl2]=ind2sub([n_theta 2*n_theta],sidx);
                clstack(k1,k2)=cl1;
                clstack(k2,k1)=cl2;
                corrstack(k1,k2)=sval;
            end          
        end

        t2=clock;
        t_cpu=t_cpu+etime(t2,t1);

    end
end


if verbose_progress
    fprintf('\n');
end

clstack_gpu=gather(g_clstack);
disp(norm(clstack(:)-clstack_gpu(:))/norm(clstack(:)))
corrstack_gpu=gather(g_corrstack);
disp(norm(corrstack(:)-corrstack_gpu(:))/norm(clstack(:)))

fprintf("t_gpu=%6.5f\n",t_gpu);
fprintf("t_cpu=%6.5f\n",t_cpu);

corrstack(corrstack~=0)=1-corrstack(corrstack~=0);
             
end
