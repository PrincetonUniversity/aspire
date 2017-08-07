function [clstack,coststack,shift_equations,shift_equations_map,clstack_mask]=...
    cryo_clmatrix_ml_ref(pf,M,KNN,noise_cov,max_itr,verbose,max_shift,shift_step,...
    ref_clmatrix,ref_shifts_2d,ref_q)
%
%   Generate common-lines matrix for the projection images' polar Fourier
%   lines given pf via maximum likelihood algorithm.
%
% The algorithm is composed of two steps: in the first step, we iteratively 
% classify the polar Fourier lines pf into M signals (clusters) via an EM 
% algorithm. That is, for each noisy line we compute soft weights indicating 
% its likelihood  of being a noisy realization of KNN signal out of the M 
% clean signals. In the second step, we detect a "hard" common line index 
% between each pair of images, by searching for two lines, one from each 
% of the two images, that attain maximum likelihood. The indices are stored
% in clstack matrix.
% Input parameters:
%
%   pf         Matrix of polar fourier rays, size n_r x n_theta x n_proj.
%   M          Number of clusters in the measurments (clean signals).
%   KNN        The number of clusters for each signal (Default = 50, 
%              must me less than M).
%   noise_cov  Noise covariance matrix, size n_r x n_r.
%   max_itr    Maximum number of iterations.
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
%   ref_q           True quaternions used to generate the projections.
% Returned variables:
%
%   clstack     Common lines matrix. (k1,k2) and (k2,k1) contain the index
%       of the common line of projections k1 and k2. (k1,k2)  contains the
%       index of the common line in the projection k1. (k2,k1) contains the
%       index of the common line in k2. 
%   coststack   The cost of the common line between projections k1
%       and k2. Since coststack is symmetric, it contain entries only above
%       the diagonal. coststack(k1,k2) measures how ''common'' is the between
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



initstate;


%% Check input parameters and set debug flags.
n_r=size(pf,1);
n_theta=size(pf,2);
n_proj=size(pf,3);
NK=n_proj;
spf=reshape(pf,n_r,n_theta*n_proj);

if (nargin<4) || (NK==-1)
    error('not enough inputs');
end

if (nargin<5)
    max_itr=0;
end

if nargin<6
    verbose=1;
end

if nargin<7
    max_shift=15; % Maximal shift between common-lines in pixels. The 
                  % shift  is from -max_shift to max_shift. 
end

if nargin<8
    shift_step=1.0; % Resolution of shift estimation in pixels.
end

max_shift=ceil(max_shift); % round up max_shift

if nargin<9
    ref_clmatrix=0;
end

if nargin<10
    ref_shifts_2d=0;
end

if nargin<11
    ref_q=0;
end

% Set flag for progress and debug messages
verbose_progress=0;
verbose_detailed_debugging=0;
verbose_plot_cl=0;
verbose_plot_shifts=0;

found_ref_clmatrix=0;
if ~isscalar(ref_clmatrix) 
    found_ref_clmatrix=1;
else
    fprintf('Reference clmatrix not found\n');
end

found_ref_shifts=0;
if ~isscalar(ref_shifts_2d)
    found_ref_shifts=1;
else
    fprintf('Reference shifts not found\n');
end

found_ref_q=0;
if ~isscalar(ref_q)
    found_ref_q=1;
else
    fprintf('Reference shifts not found\n');
end

if bitand(verbose,1)
    verbose_progress=1;
end;

if bitand(verbose,2)
        verbose_detailed_debugging=1;
end;

if bitand(verbose,4) 
    if isscalar(ref_clmatrix) 
        fprintf('Common-lines plots not available. Reference clmatrix is missing\n');
    end
    verbose_plot_cl=1;
end;

if bitand(verbose,8) 
    if isscalar(ref_clmatrix) || isscalar(ref_shifts_2d)
        fprintf('Only partial information will be plotted. Reference clmatrix or shifts are missing\n');
    end
    verbose_plot_shifts=1;
end;

verbose_plot_neighbors=0;
if bitand(verbose,16) 
    if isscalar(ref_q) 
        fprintf('Only partial information will be plotted. Reference the true projections orientation is missing\n');
    end
    verbose_plot_neighbors=1;
end;

fprintf('Verbose mode=%d\n',verbose);


    
%% Allocate variables used for estimation
msg=[];
n_proj=size(spf,2)/n_theta;
clstack=zeros(n_proj,n_proj);      % Common lines-matrix.
coststack=999*ones(n_proj,n_proj);    % Cost coefficient for each common-line.
clstack_mask=zeros(n_proj,n_proj); % Which common-lines were correctly identified.

%refdist=zeros(n_proj,n_proj); % Correlation between true common-lines.
thetadiff=zeros(n_proj,n_proj); % Angle between true and estimated common lines.

shifts_1d=zeros(n_proj,n_proj);     % Estimated 1D shift between common-lines. 

ref_shifts_1d=zeros(n_proj,n_proj); % True shift along the common-line 
    % between each pair of projections. Computed from the reference 2D
    % shifts. 
     
shift_estimation_error=zeros(n_proj,n_proj); % The difference between the 
    % estimated shift along each common line and the true shift.

% Based on the estimated common-lines, construct the equations for
% determining the 2D shift of each projection. The shift equations are
% represented using a sparse matrix, since each row in the system contains
% four non-zeros (as it involves exactly four unknowns).
% The variables below are used to construct this sparse system. The k'th
% non-zero element of the equations matrix is stored at index 
% (shift_I(k),shift_J(k)).
shift_I=zeros(4*n_proj*NK,1);  % Row index for sparse equations system.

shift_J=zeros(4*n_proj*NK,1);  % Column index for sparse equations system.

shift_eq=zeros(4*n_proj*NK,1); % The coefficients of the center estimation
    % system ordered as a single vector.
     
shift_equations_map=zeros(n_proj); % Entry (k1,k2) is the index of the 
    % euqation for the common line of projections k1 and k2. 
                               
shift_equation_idx=1;  % The equation number we are currently processing.
shift_b=zeros(n_proj*(n_proj-1)/2,1);   % Right hand side of the system.

dtheta=2*pi/n_theta;

%rk=0:n_r-1;rk=rk(:);

if verbose_progress
    fprintf('Shift estimation parameters: max_shift=%d   shift_step=%d\n',max_shift,shift_step);
end
                                                       
%% Debugging handles and variables

matched_cl=0;  % How many times the estimated common-line is close (to 
    % within a prescribed tolerance) to the true common-line.
               % Used for debugging.

if verbose_plot_cl && verbose_detailed_debugging
    h3=figure;
end
%S2CORD=0;
if verbose_plot_neighbors && found_ref_q
 %   S2CORD=R2S2(permute(q_to_rot(ref_q), [2 1 3]),n_theta); % compute S2 group of rotation matrix
end

%% Maximum likelihood based common line detector

if verbose_progress
    fprintf('Running ML...\n');
end
[S,~,Zi]=cryo_cld_ml_ws_ref(spf,n_theta,M,KNN,max_itr,noise_cov,1,max_shift,shift_step);
%[S,~,Zi]=cryo_cld_ml(spf,n_theta,M,KNN,max_itr,noise_cov,1);
%n_r=size(spf,1);
%rk=(0:(n_r-1))';
%nNN=size(Z,1);
%% Find common lines between pairs of projections
pf3=reshape(S(:,Zi(1,:)),size(spf,1),n_theta,n_proj);
rmax=size(pf3,1);
%rk=-rmax:rmax; rk=rk(:);
%n_theta=size(pf3,2);
% Prepare the shift_phases
rk2=(0:(rmax-1))';
n_shifts=ceil(2*max_shift/shift_step+1); % Number of shifts to try.
shift_phases=zeros(rmax,n_shifts);
for shiftidx=1:n_shifts
    shift=-max_shift+(shiftidx-1)*shift_step;
    shift_phases(:,shiftidx)=exp(-2*pi*sqrt(-1).*rk2.*shift./(2*rmax+1));
end
% remove the dc componenet since it is shared between the same projection image 
% and thus should not effect the detection
% for k=1:n_proj
%     pf(rmax:rmax+2,:,k)=0;
% end
corrstack=inf*ones(n_proj,n_proj);
PRECISION='single';


for k1=1:n_proj;
    
    n2=min(n_proj-k1,NK);
    subsetK2=sort(randperm(n_proj-k1)+k1);
    subsetK2=subsetK2(1:n2); % Select a subset of at most NK projections 
        % with which to search for common-lines with projection k1. 
   
    P1=pf3(:,:,k1);
    
    % Instead of processing the shifts in a loop, one shift at a time, stack
    % all shifts of a single projection into one array and compute all
    % correlations for all shifts at once.
    % The variable P1_stack stacks all shifted versions of the image k1.
    P1_stack=zeros(size(P1,1),size(P1,2)*n_shifts);
    for k=1:n_shifts
        P1_stack(:,(k-1)*n_theta+1:k*n_theta)=bsxfun(@times,P1,shift_phases(:,k));
    end

    if strcmpi(PRECISION,'single')
        g_P1_stack=gpuArray(single(P1_stack));        
    else
        g_P1_stack=gpuArray(P1_stack);
    end
    
    for k2=subsetK2;
        
        t1=clock;                       
        proj2=pf3(:,:,k2); % proj1 and proj2 are both normalized to unit norm.
        P2=proj2(1:rmax,:);
        
        if strcmpi(PRECISION,'single')
            g_P2=gpuArray(single(P2));            
        else
            g_P2=gpuArray(P2);
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
     c1=real(g_P1_stack'*g_P2);
     e1=sum(conj(g_P1_stack).*g_P1_stack);
     e1=repmat(e1',1,n_theta);
     e2=sum(conj(g_P2).*g_P2);
     e2=repmat(e2,n_theta*n_shifts,1);       
     C=e1+e2-2*c1;
    [g_sval,g_sidx]=min(C(:));
    sval=gather(g_sval);
    sidx=gather(g_sidx);
    if corrstack(k1,k2)>sval
        [cl1,sidx,cl2]=ind2sub([n_theta n_shifts n_theta],sidx);
        clstack(k1,k2)=cl1;
        clstack(k2,k1)=cl2;
        corrstack(k1,k2)=sval;
        shift=-max_shift+(sidx-1)*shift_step;
        shifts_1d(k1,k2)=shift;
        %improved_correlation=1;
    end       
%         [g_sval2,g_sidx2]=max(g_C2(:));
%         g_C1=2*real(g_P1_stack'*g_P2);
%         g_C2=2*real(g_P1_flipped_stack'*g_P2);
%         
%         [g_sval1,g_sidx1]=max(g_C1(:));
%         [g_sval2,g_sidx2]=max(g_C2(:));
       % sval1=gather(g_sval1);
       % sval2=gather(g_sval2);
       % sidx1=gather(g_sidx1);
       % sidx2=gather(g_sidx2);
        

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
        
        

        
        t2=clock;
        t=etime(t2,t1);

%         fname=sprintf('cl_4_%02d%02d',k1,k2);
%         print('-depsc',fname);


        %%%%%%%%%%%% Beginning of debug code %%%%%%%%%%%%
        % Count how many times the estimated common-line is close to within
        % max_angle of the true common-line.
        if verbose_progress
            % True line for reference
            if verbose_plot_cl
                figure(h3);
                plot(1,squeeze(corrstack(k1,k2)),'.');
                axis([0 2 0 1.1]);
            end

            % The example below assumes that the polar Fourier
            % transform was computed with n_theta=360, that is, 360
            % Fourier rays per projection. This means that C1 and C2
            % are of size 180x180. In the table below, c1 and c2 are
            % the computed common-line between proj1 and proj 2, and
            % tcl1 and tcl2 are the true common-line:
            % clstack(k1,k2,j) is always less than n_theta, so pairs of
            % common-lines that may match are (for example)
            %  cl1    cl2     tcl1  tcl2
            %   1      1        1     1
            %   1      1      181   181
            %   1    181        1   181
            %   1    181      181     1
            % These satisfy the IF statement below.
            % All the other cases
            %   1      1        1   181
            %   1      1      181     1
            %   1    181        1     1
            %   1    181      181   181
            % do not match because of orientation (one pair matches in
            % the same orientation and the other in opposite
            % orientation), and so do not statisfy the IF statement
            % below. The tables above are referred to as Table 1.

            found_matched_cl=0; % Among the NL high correlation pairs,
            % this counts the number of pairs that are close (up to
            % angle_tol) to the true common-line.

            % l1 and l2 are considered close to the true common-lines
            % tcl1 and tcl2, if the discrepancy between each line and
            % its corresponding true line is less than max_angle.
            max_angle=5/180*pi;  % 5 degrees max.
            angle_tol=2*sin(max_angle/2)+1.0e-10;

            alpha=2*pi*sqrt(-1)/(n_theta);
            PI=4*atan(1.0);

            % The estimated common-lines (l1,l2) and the true
            % common-lines (tcl1,tcl2) should be both in the same
            % orientation. For example, if the common-lines are
            % (l1,l2)=(1,1), that is the common-line between the
            % projections is in positive orientation, then the true
            % common-line is either (tcl1,tcl2)=(1,1) or
            % (tcl1,tcl2)=(180,180), that is, tcl1 and tcl2 are in the
            % same orientation. It cannot be, for example
            % (tcl1,tcl2)=(1,180). By inspecting Table 1 above we see
            % that (l1,l2) is close to (tcl1,tcl2) if
            %   a) tcll1 is close to l1 and tcll2 is close to l2. That
            %   covers the cases (l1,l2)=(1,1) (tcl1,tcl2)=(1,1), and
            %   (l1,l2)=(1,180) and (tcl1,tcl2)=(1,180). In that case,
            %   the pair (l1,l2) and the pair (tcl1,tcl2) are in the
            %   same orientation.
            %   Or
            %   b) tcl1 is close to flipped l1 and tcll2 is close to
            %   flipped l2. That covers the cases  (l1,l2)=(1,1)
            %   (tcl1,tcl2)=(180,180), and (l1,l2)=(1,180) and
            %   (tcl1,tcl2)=(180,1). In this case tcl1 is close to a
            %   flipped l1, That is, we need to flip l1 to get a
            %   match with tcl1. However, to maintain the orientation
            %   between l1 and l2 (so we stay we the same
            %   common-line), we should flip also l2, and so tcl2
            %   should match a flipped l2.
            % To conclude, for (l1,l2) to be close to (tcl1,tcl2) we
            % need that either l1 is close to tcl1 and l2 is close to
            % tcl2 (small d1s and d2s below), or, a flipped l1 is
            % close to tcl1 and a flipped l2 is close to tcl2 (small
            % d1f and d2f).
            
            if found_ref_clmatrix
                tcl1=ref_clmatrix(k1,k2);
                tcl2=ref_clmatrix(k2,k1);
                l1=clstack(k1,k2);
                l2=clstack(k2,k1);
                d1s=abs(exp(alpha*(l1-1))-exp(alpha*(tcl1-1)));
                d2s=abs(exp(alpha*(l2-1))-exp(alpha*(tcl2-1)));
                d1f=abs(exp(alpha*(l1-1)+sqrt(-1)*PI)-exp(alpha*(tcl1-1)));
                d2f=abs(exp(alpha*(l2-1)+sqrt(-1)*PI)-exp(alpha*(tcl2-1)));        
                                
                if (d1s<=angle_tol) && (d2s<=angle_tol) || ...
                        (d1f<=angle_tol) && (d2f<=angle_tol)
                    found_matched_cl=1;
                    clstack_mask(k1,k2)=1;
                    % Estimated common line is close to true common-line.
                    if verbose_plot_cl
                        hold on;
                        plot(corrstack(k1,k2),'o','MarkerSize',10,'MarkerEdgeColor','g');
                        hold off;
                    end
                else
                    % Estimated common-line is far from true common-line.
                    if verbose_plot_cl
                        hold on;
                        plot(corrstack(k1,k2),'o','MarkerSize',10,'MarkerEdgeColor','r');
                        hold off;
                    end
                end

                % Estimation error in angles
                if (tcl1<=n_theta && l1<=n_theta) ||...
                        (tcl1>n_theta && l1>n_theta)  % Same orientation for l1
                    thetadiff(k1,k2)=d1s/pi*180;
                else
                    thetadiff(k1,k2)=d1f/pi*180;
                end
                if (tcl2<=n_theta && l2<=n_theta) ||...
                        (tcl2>n_theta && l2>n_theta)  % Same orientation for l1
                    thetadiff(k2,k1)=d2s/pi*180;
                else
                    thetadiff(k2,k1)=d2f/pi*180;
                end

                if found_matched_cl
                    matched_cl=matched_cl+1;
                end
            end
        end

        %%%%%%%%%%%% End of debug code %%%%%%%%%%%%
       
         if verbose_progress                  

%              figure(3);
%              plot(triginterp(enforce_real(icfft(cryo_raynormalize(r1))),1),'b');
%              hold on;
%              plot(triginterp(enforce_real(icfft(cryo_raynormalize(r2))),1),'r');
%              hold off;
             

            if found_ref_shifts
                % Compute the true 1D shift between the common lines and
                % compare it to the  estimated shift.
                alpha=(clstack(k1,k2)-1)*dtheta;
                beta =(clstack(k2,k1)-1)*dtheta;
                dx1=ref_shifts_2d(k1,1); dy1=ref_shifts_2d(k1,2);
                dx2=ref_shifts_2d(k2,1); dy2=ref_shifts_2d(k2,2);

                % If clstack(k2,k1)==151 then beta is supposed to be exactly
                % pi. However, because of roundoff error, it might come
                % slightly less than pi, in which case the first IF will be
                % true, although it should be false. We fix this by comparing
                % beta against pi-1.0e-13 (under the reasonable assumption that
                % dtheta is larger than 1.0e-13).
                if beta<pi-1.0e-13
                    ref_shifts_1d(k1,k2)=sin(alpha)*dx1+cos(alpha)*dy1-sin(beta)*dx2-cos(beta)*dy2;
                else
                    beta=beta-pi;
                    ref_shifts_1d(k1,k2)=-sin(alpha)*dx1-cos(alpha)*dy1-sin(beta)*dx2-cos(beta)*dy2;
                end

                shift_estimation_error(k1,k2)=shifts_1d(k1,k2)-ref_shifts_1d(k1,k2);

            end

            if verbose_detailed_debugging       
                log_message('Finding common-line between projections [k1 k2]=[ %d %d ]',k1,k2);
                log_message('\t Common lines:');
                log_message('\t \t clstack= [ %3d %3d ]',clstack(k1,k2),clstack(k2,k1)); 

                if found_ref_clmatrix
                    log_message('\t \t ref    = [ %3d %3d ]',ref_clmatrix(k1,k2),ref_clmatrix(k2,k1));
                    log_message('\t \t dtheta1= %5.2f (degrees)',thetadiff(k1,k2));
                    log_message('\t \t dtheta2= %5.2f (degrees)',thetadiff(k2,k1));
                    if found_matched_cl                       
                        log_message('\t \t status = MATCHED (less than %5.2f degrees)',max_angle/pi*180);
                    else
                        log_message('\t \t status = NOT MATCHED (more than %5.2f degrees)',max_angle/pi*180');
                    end
                end

                log_message('\t Correlation:');
                log_message('\t \t est    =  %9.6f',corrstack(k1,k2));

                if found_ref_clmatrix && found_ref_shifts
              %      log_message('\t \t true   =  %9.6f',refcorr(k1,k2));
                end

                log_message('\t Shifts:');
                log_message('\t \t est  shift= %6.4f',shifts_1d(k1,k2));

                if found_ref_shifts
                    log_message('\t \t True shift= %6.4f',ref_shifts_1d(k1,k2));
                    log_message('\t \t shift err = %3.2e',shift_estimation_error(k1,k2));
                end
                log_message(' ');
            end
        end

        % Create a shift equation for the projections pair (k1,k2).
        idx=4*(shift_equation_idx-1)+1:4*shift_equation_idx;
        shift_alpha=(clstack(k1,k2)-1)*dtheta;  % Angle of common ray in projection 1.
        shift_beta= (clstack(k2,k1)-1)*dtheta;  % Angle of common ray in projection 2.
        shift_I(idx)=shift_equation_idx; % Row index to construct the sparse equations.
        shift_J(idx)=[2*k1-1 2*k1 2*k2-1 2*k2]; % Columns of the shift variables that correspond to the current pair (k1,k2).
        shift_b(shift_equation_idx)=shifts_1d(k1,k2); % Right hand side of the current equation.

        % Compute the coefficients of the current equation.
        if shift_beta<pi       
            shift_eq(idx)=[sin(shift_alpha) cos(shift_alpha) -sin(shift_beta) -cos(shift_beta)];
        else
            shift_beta=shift_beta-pi; % In the derivation we assume that all angles are less 
                                      % than PI where angles larger than PI are assigned 
                                      % nigative orientation.
            shift_eq(idx)=[-sin(shift_alpha) -cos(shift_alpha) -sin(shift_beta) -cos(shift_beta)];
        end
    
        shift_equations_map(k1,k2)=shift_equation_idx;  % For each pair (k1,k2), store the index of its equation.
        shift_equation_idx=shift_equation_idx+1;
        
        
        if verbose_progress
            bs=char(repmat(8,1,numel(msg)));
            fprintf('%s',bs);
            msg=sprintf('k1=%3d/%3d  k2=%3d/%3d  t=%7.5f',k1,n_proj,k2,n_proj,t);
            fprintf('%s',msg);
        end
    end;
end;














% 
% for k1=1:n_proj;
%     t1=clock; 
%     n2=min(n_proj-k1,NK);
%     subsetK2=sort(randperm(n_proj-k1)+k1);
%     subsetK2=subsetK2(1:n2); % Select a subset of at most NK projections 
%         % with which to search for common-lines with projection k1.
%     P1_CL=Zi(:,(k1-1)*n_theta+(1:n_theta));
%     Z1=Z(:,(k1-1)*n_theta+(1:n_theta)); 
%     for k2=subsetK2;
%          P2_CL=Zi(:,(k2-1)*n_theta+(1:n_theta));
%          Z2=Z(:,(k2-1)*n_theta+(1:n_theta));
%          [~,cl1,cl2]=intersect(P1_CL(:),P2_CL(:)); % Search pair of lines (from projection
%         %images P1 and P2) sharing the same cluster.
%        
%        if ~isempty(cl1)
%         %   compute the maximum likelihood 
%             Z_cl1=Z1(cl1);
%             Z_cl2=Z2(cl2);
%             Z_cl=Z_cl1.*Z_cl2;
%             [cl_cost,i]=max(Z_cl);
%             cl_shift_1d=0;%SHIFT(i);
%             cl1=ceil(cl1(i)/nNN);
%             cl2=ceil(cl2(i)/nNN);
%        else
%             T=floor(max_shift/shift_step);
%             shifts=(-T:T)*shift_step;
%             
%             cl_cost=inf*ones(1,length(shifts));
%             cl_idx=zeros(1,length(shifts));
%             for shiftidx=1:length(shifts) 
%                 shift_phases=exp(2*pi*sqrt(-1).*rk.*shifts(shiftidx)./(2*n_r+1)); 
%                 shift_phases=repmat(shift_phases,1,n_theta);
%                 P1_shifted=pf3(:,:,k1).*shift_phases;  
%                 P2=pf3(:,:,k2); 
%                 c1=real(P1_shifted'*P2);
%                 e1=sum(conj(P1_shifted).*P1_shifted);
%                 e1=repmat(e1',1,n_theta);
%                 e2=sum(conj(P2).*P2);
%                 e2=repmat(e2,n_theta,1);       
%                 % compute the distance cost between P1 and P2
%                 C=e1+e2-2*c1;
%                 % The best match is the minimum 
%                 [LL_cl,sidx]=min(C(:));
%                 cl_cost(shiftidx)=LL_cl;
%                 cl_idx(shiftidx)=sidx;
%                 %if LL_cl<min(cl_cost)
%                     %cl_cost(shiftidx)=LL_cl;
%                     %[cl1,cl2]=ind2sub([n_theta n_theta],sidx);
%                     %improved_correlation=1;
%                     %cl_shift_1d=shifts(shiftidx);
%                 %end
%             end
%             [~,clx]=min(cl_cost);
%             cl_shift_1d=shifts(clx);
%             sidx=cl_idx(clx);
%             [cl1,cl2]=ind2sub([n_theta n_theta],sidx);
%        end
%         
%         clstack(k1,k2)=cl1;
%         clstack(k2,k1)=cl2;
%         coststack(k1,k2)=cl_cost;
%         shifts_1d(k1,k2)=cl_shift_1d;
% 
%         %%%%%%%%%%%% Beginning of debug code %%%%%%%%%%%%
%         % Count how many times the estimated common-line is close to within
%         % max_angle of the true common-line.
%         if verbose_detailed_debugging
%             % True line for reference
%             if verbose_plot_cl
%                 figure(h3);
%                 plot(1,squeeze(coststack(k1,k2)),'.');
%                 axis([0 2 0 1.1]);
%             end
% 
%             % The example below assumes that the polar Fourier
%             % transform was computed with n_theta=360, that is, 360
%             % Fourier rays per projection. This means that C
%             % is of size 360x360. In the table below, cl1 and cl2 are
%             % the computed common-line between proj1 and proj 2, and
%             % tcl1 and tcl2 are the true common-line:
%             % clstack(k1,k2,j) is always less than n_theta, so pairs of
%             % common-lines that may match are (for example)
%             %  cl1    cl2     tcl1  tcl2
%             %   1      1        1     1
%             %   1      1      181   181
%             %   1    181        1   181
%             %   1    181      181     1
%             % These satisfy the IF statement below.
%             % All the other cases
%             %   1      1        1   181
%             %   1      1      181     1
%             %   1    181        1     1
%             %   1    181      181   181
%             % do not match because of orientation (one pair matches in
%             % the same orientation and the other in opposite
%             % orientation), and so do not statisfy the IF statement
%             % below. The tables above are referred to as Table 1.
% 
%             found_matched_cl=0; % Among the NL high correlation pairs,
%             % this counts the number of pairs that are close (up to
%             % angle_tol) to the true common-line.
% 
%             % l1 and l2 are considered close to the true common-lines
%             % tcl1 and tcl2, if the discrepancy between each line and
%             % its corresponding true line is less than max_angle.
%             max_angle=5/180*pi;  % 5 degrees max.
%             angle_tol=2*sin(max_angle/2)+1.0e-10;
% 
%             alpha=2*pi*sqrt(-1)/n_theta;
%             PI=4*atan(1.0);
% 
%             % The estimated common-lines (l1,l2) and the true
%             % common-lines (tcl1,tcl2) should be both in the same
%             % orientation. For example, if the common-lines are
%             % (l1,l2)=(1,1), that is the common-line between the
%             % projections is in positive orientation, then the true
%             % common-line is either (tcl1,tcl2)=(1,1) or
%             % (tcl1,tcl2)=(180,180), that is, tcl1 and tcl2 are in the
%             % same orientation. It cannot be, for example
%             % (tcl1,tcl2)=(1,180). By inspecting Table 1 above we see
%             % that (l1,l2) is close to (tcl1,tcl2) if
%             %   a) tcll1 is close to l1 and tcll2 is close to l2. That
%             %   covers the cases (l1,l2)=(1,1) (tcl1,tcl2)=(1,1), and
%             %   (l1,l2)=(1,180) and (tcl1,tcl2)=(1,180). In that case,
%             %   the pair (l1,l2) and the pair (tcl1,tcl2) are in the
%             %   same orientation.
%             %   Or
%             %   b) tcl1 is close to flipped l1 and tcll2 is close to
%             %   flipped l2. That covers the cases  (l1,l2)=(1,1)
%             %   (tcl1,tcl2)=(180,180), and (l1,l2)=(1,180) and
%             %   (tcl1,tcl2)=(180,1). In this case tcl1 is close to a
%             %   flipped l1, That is, we need to flip l1 to get a
%             %   match with tcl1. However, to maintain the orientation
%             %   between l1 and l2 (so we stay we the same
%             %   common-line), we should flip also l2, and so tcl2
%             %   should match a flipped l2.
%             % To conclude, for (l1,l2) to be close to (tcl1,tcl2) we
%             % need that either l1 is close to tcl1 and l2 is close to
%             % tcl2 (small d1s and d2s below), or, a flipped l1 is
%             % close to tcl1 and a flipped l2 is close to tcl2 (small
%             % d1f and d2f).
%             
%             if found_ref_clmatrix
%                 tcl1=ref_clmatrix(k1,k2);
%                 tcl2=ref_clmatrix(k2,k1);
%                 l1=clstack(k1,k2);
%                 l2=clstack(k2,k1);
%                 d1s=abs(exp(alpha*(l1-1))-exp(alpha*(tcl1-1)));
%                 d2s=abs(exp(alpha*(l2-1))-exp(alpha*(tcl2-1)));
%                 d1f=abs(exp(alpha*(l1-1)+sqrt(-1)*PI)-exp(alpha*(tcl1-1)));
%                 d2f=abs(exp(alpha*(l2-1)+sqrt(-1)*PI)-exp(alpha*(tcl2-1)));
% 
%                 if (d1s<=angle_tol) && (d2s<=angle_tol) || ...
%                         (d1f<=angle_tol) && (d2f<=angle_tol)
%                     found_matched_cl=1;
% 
%                     % Estimated common line is close to true common-line.
%                     if verbose_plot_cl
%                         %hold on;
%                         %plot(coststack(k1,k2),'o','MarkerSize',10,'MarkerEdgeColor','g');
%                         %hold off;
%                     end
%                 else
%                     % Estimated common-line is far from true common-line.
%                     if verbose_plot_cl
%                         %hold on;
%                         %plot(coststack(k1,k2),'o','MarkerSize',10,'MarkerEdgeColor','r');
%                         %hold off;
%                     end
%                 end
% 
%                 % Estimation error in angles
%                 if (tcl1<=n_theta && l1<=n_theta) ||...
%                         (tcl1>n_theta && l1>n_theta)  % Same orientation for l1
%                     thetadiff(k1,k2)=d1s/pi*180;
%                 else
%                     thetadiff(k1,k2)=d1f/pi*180;
%                 end
%                 if (tcl2<=n_theta && l2<=n_theta) ||...
%                         (tcl2>n_theta && l2>n_theta)  % Same orientation for l1
%                     thetadiff(k2,k1)=d2s/pi*180;
%                 else
%                     thetadiff(k2,k1)=d2f/pi*180;
%                 end
% 
%                 if found_matched_cl
%                     matched_cl=matched_cl+1;
%                 end
%             end
% 
% %             if verbose_detailed_debugging && found_ref_clmatrix && found_ref_shifts
% %                 % Compute and store the distance between the true
% %                 % common-lines.
% %                 l1=ref_clmatrix(k1,k2);
% %                 l2=ref_clmatrix(k2,k1);
% % 
% %                 r1=P1(:,l1);
% %                 r2=P2(:,l2);
% %                 
% % 
% %                 alpha=(l1-1)*dtheta;
% %                 beta =(l2-1)*dtheta;
% %                 dx1=ref_shifts_2d(k1,1); dy1=ref_shifts_2d(k1,2);
% %                 dx2=ref_shifts_2d(k2,1); dy2=ref_shifts_2d(k2,2);
% %                 % Shift by the exact amount:
% %                 ds=sin(alpha)*dx1+cos(alpha)*dy1-sin(beta)*dx2-cos(beta)*dy2;
% %                 phi=exp(-2*pi*sqrt(-1).*rk.*ds./(2*n_r+1));
% %                 r1=r1.*phi;
% %                 refdist(k1,k2)=enforce_real((r1-r2)'*(r1-r2));
% %             end
% 
%             %%%%%%%%%%%% Beginning of debug code %%%%%%%%%%%%            
%             if verbose_plot_shifts && improved_correlation
%             % If the current shift produces a better correlation, then plot
%             % the two appropriately shifted estimated common-lines. The
%             % figure also displays the true and estimated common lines, and
%             % the estimated correlation value.
%             
% %                 fac=max(1/shift_step,1);
% %                 xx=-rmax:1/fac:rmax;
% %                 
% %                 if clstack(k2,k1)<=n_theta                               
% %                     proj1_shifted=[P1; zeros(1,n_theta) ; conj(flipud(P1))];
% %                     p1=enforce_real(icfft(proj1_shifted(:,clstack(k1,k2))));                                                            
% %                     p2=enforce_real(icfft(P2(:,clstack(k2,k1))));
% %                 else
% %                     proj1_shifted_flipped=[P1_shifted_flipped; zeros(1,n_theta) ; conj(flipud(P1_shifted_flipped))];
% %                     p1=enforce_real(icfft(proj1_shifted_flipped(:,clstack(k1,k2))));                   
% %                     p2=enforce_real(icfft(P2(:,clstack(k2,k1)-n_theta)));
% %                 end
%                 
%                 %v1=triginterp(p1,fac);
%                 %v2=triginterp(p2,fac);
%                 
%                 %figure(h2);
%                 %plot(xx,real(v1));
%                 %hold on;                
%                 %plot(xx,real(v2),'r');
%                 %px=get(gca,'Xlim');
%                 %py=get(gca,'Ylim');
%                                 
%                 if found_ref_clmatrix
%                     l1str=sprintf('%3d',ref_clmatrix(k1,k2));
%                     l2str=sprintf('%3d',ref_clmatrix(k2,k1));
%                 else
%                     l1str='N/A';
%                     l2str='N/A';
%                 end
%                 
%                 str=sprintf(strcat(' k1=%d  k2=%d \n shift = %4.2f \n',...
%                     ' est cl   = [ %3d  %3d ] \n',...
%                     ' ref cl   = [ %s  %s ] \n',...
%                     ' est corr = %7.5f'),...
%                     k1,k2,shift,clstack(k1,k2),clstack(k2,k1),...
%                     l1str,l2str,...
%                     coststack(k1,k2));
%                 
%                 if found_ref_clmatrix && found_ref_shifts
%                     str=strcat(str,...
%                         sprintf('\n ref corr = %7.5f',refdist(k1,k2)));
%                 end
% 
%                 text(px(2)*0.4,py(2)*0.7,str,'EdgeColor','k')
%                 hold off;                                
%             end;
%             %%%%%%%%%%%% End of debug code %%%%%%%%%%%%
%         end
% 
%         t2=clock;
%         t=etime(t2,t1);
% 
%         %%%%%%%%%%%% End of debug code %%%%%%%%%%%%
%        
%          if verbose_detailed_debugging                        
% 
%             if found_ref_shifts
%                 % Compute the true 1D shift between the common lines and
%                 % compare it to the  estimated shift.
%                 alpha=(clstack(k1,k2)-1)*dtheta;
%                 beta =(clstack(k2,k1)-1)*dtheta;
%                 dx1=ref_shifts_2d(k1,1); dy1=ref_shifts_2d(k1,2);
%                 dx2=ref_shifts_2d(k2,1); dy2=ref_shifts_2d(k2,2);
% 
%                 % If clstack(k2,k1)==151 then beta is supposed to be exactly
%                 % pi. However, because of roundoff error, it might come
%                 % slightly less than pi, in which case the first IF will be
%                 % true, although it should be false. We fix this by comparing
%                 % beta against pi-1.0e-13 (under the reasonable assumption that
%                 % dtheta is larger than 1.0e-13).
%                 if beta<pi-1.0e-13
%                     ref_shifts_1d(k1,k2)=sin(alpha)*dx1+cos(alpha)*dy1-sin(beta)*dx2-cos(beta)*dy2;
%                 else
%                     beta=beta-pi;
%                     ref_shifts_1d(k1,k2)=-sin(alpha)*dx1-cos(alpha)*dy1-sin(beta)*dx2-cos(beta)*dy2;
%                 end
% 
%                 shift_estimation_error(k1,k2)=shifts_1d(k1,k2)-ref_shifts_1d(k1,k2);
% 
%             end
% 
%             if found_ref_clmatrix
% %               
%                 if found_matched_cl
%                     clstack_mask(k1,k2)=1;
%                 end
%             end
% 
%         end
% 
%         % Create a shift equation for the projections pair (k1,k2).
%         idx=4*(shift_equation_idx-1)+1:4*shift_equation_idx;
%         shift_alpha=(clstack(k1,k2)-1)*dtheta;  % Angle of common ray in projection 1.
%         shift_beta= (clstack(k2,k1)-1)*dtheta;  % Angle of common ray in projection 2.
%         shift_I(idx)=shift_equation_idx; % Row index to construct the sparse equations.
%         shift_J(idx)=[2*k1-1 2*k1 2*k2-1 2*k2]; % Columns of the shift variables that correspond to the current pair (k1,k2).
%         shift_b(shift_equation_idx)=shifts_1d(k1,k2); % Right hand side of the current equation.
% 
%         % Compute the coefficients of the current equation.
%         if shift_beta<pi       
%             shift_eq(idx)=[sin(shift_alpha) cos(shift_alpha) -sin(shift_beta) -cos(shift_beta)];
%         else
%             shift_beta=shift_beta-pi; % In the derivation we assume that all angles are less 
%                                       % than PI where angles larger than PI are assigned 
%                                       % nigative orientation.
%             shift_eq(idx)=[-sin(shift_alpha) -cos(shift_alpha) -sin(shift_beta) -cos(shift_beta)];
%         end
%     
%         shift_equations_map(k1,k2)=shift_equation_idx;  % For each pair (k1,k2), store the index of its equation.
%         shift_equation_idx=shift_equation_idx+1;
%         
%         
%         if verbose_progress
%             bs=char(repmat(8,1,numel(msg)));
%             fprintf('%s',bs);
%             msg=sprintf('k1=%3d/%3d  k2=%3d/%3d  t=%7.5f\n',k1,n_proj,k2,n_proj,t);
%             fprintf('%s',msg);
%         end
%    
%     end;
% end;


if verbose_progress
    fprintf('\n');
end

if verbose_detailed_debugging && found_ref_clmatrix
   fprintf('Matched common-lines=%d/%d  (%3.1f%%)\n',...
       matched_cl,n_proj*(n_proj-1)/2,100*matched_cl/(((n_proj*(n_proj-1))/2)));
end

             
% Construct least-squares for the two-dimensioal shifts.
shift_equation_idx=shift_equation_idx-1;
shift_equations=sparse(shift_I(1:4*shift_equation_idx),...
    shift_J(1:4*shift_equation_idx),shift_eq(1:4*shift_equation_idx),...
    shift_equation_idx,2*n_proj);

shift_equations=[shift_equations shift_b(1:shift_equation_idx)];


if verbose_detailed_debugging
    % XXX Check that the shift estimation improves with n_theta.
    % XXX Does it improve with n_proj?
    if n_proj<=100
        [~,~,V]=svd(full(shift_equations(:,1:end-1)));
        %s=diag(S); % take the spectrum
        %fprintf('Singular values of the shift system of equations:');
        %fprintf('%d  \n',fliplr(s.'));

        % Check that the difference between the true shifts and the estimated ones
        % is in the null space of the equations.
        est_shifts=shift_equations(:,1:end-1)\shift_equations(:,end);
        est_shifts=transpose(reshape(est_shifts,2,n_proj));
        
        if found_ref_shifts
            s1=reshape(ref_shifts_2d.',2*n_proj,1);
            s2=reshape(est_shifts.',2*n_proj,1);
            V=V(:,1:end-3); % Null space of shift_equations.
            % Compute the difference between the true and estimated shifts in
            % the subspace that is orthogonal to the null space of
            % shift_equations.

            if norm(V.'*s1)>1.0e-12
                fprintf('Difference between true and estimated shifts: %8.5e\n',...
                    (norm(V.'*(s1-s2))/norm(V.'*s1)));
            else
                fprintf('norm(V.''*s1) = %7.5e\n',norm(V.'*s1));
            end
        end
    else
        fprintf('Not computing SVD of shifts matrix -- matrix is too big\n');
    end;
end



