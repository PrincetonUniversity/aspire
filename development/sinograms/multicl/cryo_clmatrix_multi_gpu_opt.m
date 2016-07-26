function [cl_multi_stack,corr_multi_stack,shifts_1d_multi]=...
    cryo_clmatrix_multi_gpu_opt(pf,NK,max_shift,shift_step,max_angle, Ncl,multisuppression)

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
%   verbose  Bitmask of debugging level.Bits:
%           0   silent
%           1   One line progress message (not written to log) (Default)
%           2   Print detailed debug messages
%           4   Draw common-line debugging plots
%           8   Draw shift estimation debugging plots
%   max_shift       Maximal 1D shift (in pixels)  to search between
%       common-lines. Default: 15.
%   shift_step      Resolution of shift estimation in pixels. Note that
%        shift_step can be any positive real number. Default: 1.
%   ref_clmatrix    True common-lines matrix (for debugging).
%   ref_shifts_2d   True 2D shifts between projections (for debugging).
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
%   clstack_mask  If ref_clmatrix is given, this array discribes which
%       common-lines were identified correcly. It is of size
%       n_projXn_projs, where entry (k1,k2) is 1 if the common-line between
%       projections k1 and k2 was correctly identified, and 0 otherwise.
%       This matrix will be non-zero only if bit 2 of verbose it set.
%
% Yoel Shkolnisky, May 2013.
%
% Revision:
% 16/05/2013 16:54 Yoel Shkolnisky
%         Optimize code:
%             1. Compute the phases required for all shifts only once and
%             not inside the shifts loop.
%             2. Add a PRECISION variable that determines the precision
%             (single/double) in which calculations are done in the GPU.
%             3. To minimize communication between the CPU and GPU, search
%             for the best correlation inside the GPU. and transfer to the
%             CPU only the value of the best correlation and the index of
%             the best common line and its shift.
%             4. The debug code inside the shifts loop (look for "Beginning
%             of debug code") does not work any more, since clstack is not
%             being updated inside the loop any more. This debug code
%             should be moved outside the shifts loop (in future versions).
%
%  16/05/2013 20:33 Yoel Shkolnisky
%          Further optimization: Instead of processing the shifts in a
%          loop, one shift at a time, stack all shifts of a single
%          projection into one array and compute all correlations for all
%          shifts at once.
%
%  16/05/2013 21:00 Yoel Shkolnisky
%          Code cleanup. Print warning when using single precision GPU
%          calulations. Speedup of this function vs cryo_clmatrix is ~4.5
%          for double precision and ~14 for single precision.



PRECISION='single';

%log_message('GPU PRECISION=%s',PRECISION);
% if strcmpi(PRECISION,'single')
%     log_message('Using single precision GPU calculations!!!');
% end

%initstate;%why this is her?
msg=[];

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

% XXX The PCA should not be done here. Move outside the common lines
% XXX matrix.
% Project all common lines on the first 15 principal components.
% % % nPCs=30;
% % % spf=reshape(pf,size(pf,1),size(pf,2)*size(pf,3));
% % % 
% % % m=mean(spf,2);
% % % spf1=zeros(size(spf));
% % % for k=1:size(spf,2);
% % %     spf1(:,k)=spf(:,k)-m;
% % % end
% % % spf=spf1;
% % % 
% % % [U,S]=svd(spf);
% % % U=U(:,1:nPCs);
% % % pf_sav=pf;
% % % spf=U*U'*spf;
% % % pf=reshape(spf,size(pf,1),size(pf,2),size(pf,3));
% % % 
% % % log_message('Number of PCs = %d',nPCs);

n_theta=size(pf,2);
n_proj=size(pf,3);

%% Check input parameters and set debug flags.
if (nargin<2) || (NK==-1)
    NK=n_proj; % Number of common-line pairs to compute for each projection
end


if nargin<3
    max_shift=15; % Maximal shift between common-lines in pixels. The
    % shift  is from -max_shift to max_shift.
end

if nargin<4
    shift_step=1.0; % Resolution of shift estimation in pixels.
end
n_shifts=2*max_shift/shift_step+1; % Number of shifts to try.



%new for multi
%%


if nargin<5
    max_angle=5;  % 5 degrees max.
end

if nargin<6
    Ncl=10;
end

if nargin<7
    multisuppression=1; %we suppress maximum in all shifts
end

%%

% Set flag for progress and debug messages

verbose_plot_shifts=0;


found_ref_clmatrix=0;


%%
%new for multi
peakjump=round( max_angle/ (180 / n_theta));
%%


%%
%new for multi
cl_multi_stack=cell(Ncl,1);
cl_multi_stack(:)={zeros(n_proj,n_proj,PRECISION)}; %Common lines-matrix with 5 common lines for each pair of projections taken from different peaks in the hisogram.


corr_multi_stack=cell(Ncl,1);
corr_multi_stack(:)={zeros(n_proj)}; % Correlation coefficient for each common-line (5 cl for each pair of projections [not from same peak]).


%% Allocate variables used for shift estimation
%%
%new for multi
shifts_1d_multi=cell(Ncl,1);
shifts_1d_multi(:)={zeros(n_proj,n_proj,PRECISION)};




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

rk2=rk(1:rmax);
% Prepare the shift_phases
shift_phases=zeros(rmax,n_shifts);
for shiftidx=1:n_shifts
    shift=-max_shift+(shiftidx-1)*shift_step;
    shift_phases(:,shiftidx)=exp(-2*pi*sqrt(-1).*rk2.*shift./(2*rmax+1));
end


for k1=1:n_proj;
    
    n2=min(n_proj-k1,NK);
    subsetK2=sort(randperm(n_proj-k1)+k1);
    subsetK2=subsetK2(1:n2); % Select a subset of at most NK projections
    % with which to search for common-lines with projection k1.
    
    proj1=pf3(:,:,k1);
    P1=proj1(1:rmax,:);  % Take half ray plus the DC
    P1_flipped=conj(P1);
    
    % Instead of processing the shifts in a loop, one shift at a time, stack
    % all shifts of a single projection into one array and compute all
    % correlations for all shifts at once.
    % The variable P1_stack stacks all shifted versions of the image k1.
    P1_stack=zeros(size(P1,1),size(P1,2)*n_shifts);
    P1_flipped_stack=zeros(size(P1,1),size(P1,2)*n_shifts);
    for k=1:n_shifts
        P1_stack(:,(k-1)*n_theta+1:k*n_theta)=bsxfun(@times,P1,shift_phases(:,k));
        P1_flipped_stack(:,(k-1)*n_theta+1:k*n_theta)=bsxfun(@times,P1_flipped,shift_phases(:,k));
    end
    
    if strcmpi(PRECISION,'single')
        g_P1_stack=gpuArray(single(P1_stack));
        g_P1_flipped_stack=gpuArray(single(P1_flipped_stack));
    else
        g_P1_stack=gpuArray(P1_stack);
        g_P1_flipped_stack=gpuArray(P1_flipped_stack);
    end
    
    
    % Make sure the DC component is zero. This is assumed  below in
    % computing correlations.
    if norm(proj1(rmax+1,:))>1.0e-13
        error('DC component of projection is not zero');
    end
    
    for k2=subsetK2;
        
        
        proj2=pf3(:,:,k2); % proj1 and proj2 are both normalized to unit norm.
        P2=proj2(1:rmax,:);
        
        if strcmpi(PRECISION,'single')
            g_P2=gpuArray(single(P2));
        else
            g_P2=gpuArray(P2);
        end
        
        if norm(proj2(rmax+1,:))>1.0e-13
            error('DC component of projection is not zero');
        end
        
        
        %%%%%%%%%%%% Beginning of debug code %%%%%%%%%%%%
        if verbose_plot_shifts && found_ref_clmatrix
            % Plot the signals that correspond to the true common-line.
            % This allows to appreciate visually that this is indeed the
            % common line, as well as the shift between the signal. Note
            % that the plotted signal are not filtered.
            xx=-rmax:1:rmax;
            if ref_clmatrix(k1,k2)<=n_theta
                pf1=pf(:,ref_clmatrix(k1,k2),k1);
            else
                pf1=pf(:,ref_clmatrix(k1,k2)-n_theta,k1);
                pf1=flipud(pf1);
            end
            
            p1=enforce_real(cfft(cryo_raynormalize(pf1)));
            v1=triginterp(p1,1);
            
            if ref_clmatrix(k2,k1)<=n_theta
                pf2=pf(:,ref_clmatrix(k2,k1),k2);
            else
                pf2=pf(:,ref_clmatrix(k2,k1)-n_theta,k2);
                pf2=flipud(pf2);
            end
            
            p2=enforce_real(cfft(cryo_raynormalize(pf2)));
            v2=triginterp(p2,1);
            
            figure(h1);
            plot(xx,real(v1),'-b');
            hold on;
            plot(xx,real(v2),'-r');
            hold off;
            % We always measure by how much we need to shift the blue
            % signal (v1). If we need to shift it to the right then shift is
            % negative and dx is positive.
            legend('Shifted signal','Fixed Signal')
            
        end
        %%%%%%%%%%%% End of debug code %%%%%%%%%%%%
        
        
        %Optimized GPU version:
        g_C1=2*real(g_P1_stack'*g_P2);
        g_C2=2*real(g_P1_flipped_stack'*g_P2);
        
        %         [g_sval1,g_sidx1]=max(g_C1(:));
        %         [g_sval2,g_sidx2]=max(g_C2(:));
        %         sval1=gather(g_sval1);
        %         sval2=gather(g_sval2);
        %         sidx1=gather(g_sidx1);
        %         sidx2=gather(g_sidx2);
        
        
        % Optimized CPU version. Uncommet to use on CPU:
        %         C1=2*real(P1_stack'*P2);
        %         C2=2*real(P1_flipped_stack'*P2);
        %
        %         [sval1,sidx1]=max(C1(:));
        %         [sval2,sidx2]=max(C2(:));
        
        %         if sval1>sval2
        %             [cl1,sidx,cl2]=ind2sub([n_theta n_shifts n_theta],sidx1);
        %             clstack(k1,k2)=cl1;
        %             clstack(k2,k1)=cl2;
        %             corrstack(k1,k2)=sval1;
        %             shift=-max_shift+(sidx-1)*shift_step;
        %             shifts_1d(k1,k2)=shift;
        %             improved_correlation=1;
        %         else
        %             [cl1,sidx,cl2]=ind2sub([n_theta n_shifts n_theta],sidx2);
        %             clstack(k1,k2)=cl1;
        %             clstack(k2,k1)=cl2+n_theta;
        %             corrstack(k1,k2)=sval2;
        %             shift=-max_shift+(sidx-1)*shift_step;
        %             shifts_1d(k1,k2)=shift;
        %             improved_correlation=1;
        %         end
        
        
        %%%%%%%%%%%%%%%%%%%%% new for multi BEGIN
        %%
        
        
        g_C1peak=g_C1;
        g_C2peak=g_C2;
        
        for cl=1:Ncl
            
            
            [g_sval1,g_sidx1]=max(g_C1peak(:));
            [g_sval2,g_sidx2]=max(g_C2peak(:));
            %    sval1=gather(g_sval1);
            %    sval2=gather(g_sval2);
            %    sidx1=gather(g_sidx1);
            %    sidx2=gather(g_sidx2);
            
            G = gather([g_sval1, g_sval2, g_sidx1, g_sidx2]);
            [sval1, sval2, sidx1, sidx2]=deal(G(1),G(2),G(3),G(4));
            if sval1>sval2
                [cl1,sidx,cl2]=ind2sub([n_theta n_shifts n_theta],sidx1);
                cl_multi_stack{cl}(k1,k2)=cl1;
                cl_multi_stack{cl}(k2,k1)=cl2;
                corr_multi_stack{cl}(k1,k2)=sval1;
                shift=-max_shift+(sidx-1)*shift_step;
                shifts_1d_multi{cl}(k1,k2)=shift;
                adgpeakjump=peakjump;
                while (cl1<=adgpeakjump || cl2<=adgpeakjump || cl1>n_theta-adgpeakjump || cl2>n_theta-adgpeakjump)
                    adgpeakjump=floor(adgpeakjump/2);
                end
                if (multisuppression==1)
                    
                    suppressioncolmn=ones(n_theta*n_shifts,adgpeakjump*2+1);
                    for sidxx=1:n_shifts
                        suppressioncolmn((sidxx-1)*n_theta+cl1-adgpeakjump:(sidxx-1)*n_theta+cl1+adgpeakjump,:)=0;
                        
                    end
                    g_C1peak(:, cl2-adgpeakjump:cl2+adgpeakjump)=...
                        g_C1peak(:, cl2-adgpeakjump:cl2+adgpeakjump).*suppressioncolmn;
                    
                    
                    
                else
                    g_C1peak((sidx-1)*n_theta+cl1-adgpeakjump:(sidx-1)*n_theta+cl1+adgpeakjump,...
                        cl2-adgpeakjump:cl2+adgpeakjump)=0;
                end
            else
                [cl1,sidx,cl2]=ind2sub([n_theta n_shifts n_theta],sidx2);
                cl_multi_stack{cl}(k1,k2)=cl1;
                cl_multi_stack{cl}(k2,k1)=cl2+n_theta;
                corr_multi_stack{cl}(k1,k2)=sval2;
                shift=-max_shift+(sidx-1)*shift_step;
                shifts_1d_multi{cl}(k1,k2)=shift;
                adgpeakjump=peakjump;
                while (cl1<=adgpeakjump || cl2<=adgpeakjump || cl1>n_theta-adgpeakjump || cl2>n_theta-adgpeakjump)
                    adgpeakjump=floor(adgpeakjump/2);
                end
                if (multisuppression==1)
                    
                    
                    suppressioncolmn=ones(n_theta*n_shifts,adgpeakjump*2+1);
                    for sidxx=1:n_shifts
                        suppressioncolmn((sidxx-1)*n_theta+cl1-adgpeakjump:(sidxx-1)*n_theta+cl1+adgpeakjump,:)=0;
                        
                    end
                    g_C2peak(:, cl2-adgpeakjump:cl2+adgpeakjump)=...
                        g_C2peak(:, cl2-adgpeakjump:cl2+adgpeakjump).*suppressioncolmn;
                    
                    
                else
                    g_C2peak((sidx-1)*n_theta+cl1-adgpeakjump:(sidx-1)*n_theta+cl1+adgpeakjump,...
                        cl2-adgpeakjump:cl2+adgpeakjump)=0;
                end
                
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%% new for multi END
%%


for nclind=1:Ncl
    tmp=corr_multi_stack{nclind};
    tmp(tmp~=0)=1-tmp(tmp~=0);
    corr_multi_stack{nclind}=tmp;
end


        %         fname=sprintf('cl_4_%02d%02d',k1,k2);
        %         print('-depsc',fname);
        
        
       
       

function y=enforce_real(x)
% The polar Fourier transform of each projection is computed to single
% precision. The correlations should therefore be real to single precision.
err=norm(imag(x(:)))/norm(x(:));
if err>1.0e-7
    warning('GCAR:imaginaryComponents','Imaginary components err=%d',err);
end
y=real(x);
