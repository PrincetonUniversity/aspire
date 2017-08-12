function [est_shifts,shift_equations]=cryo_estimate_shifts(pf,rotations,...
    max_shift,shift_step,memoryfactor,shifts_2d_ref,verbose)
%
% CRYO_ESTIMATE_SHIFTS  Esimate 2D shifts in projections.
%
% est_shifts=cryo_estimate_shifts(pf,rotations,max_shift,shift_step,shifts_2d_ref)
%   Estimate 2D shifts in projections using the estimated rotations.
%   The function computes the common lines from the estimated rotations,
%   and then, for each common line, estimates the 1D shift between its two
%   Fourier rays (one in image i and one in image j). Using the common
%   lines and the 1D shifts, the function solves the least-squares
%   equations for the 2D shifts.
%
%   The function processes the (Fourier transformed) projections exactly as
%   cryo_clmatrix.
%
% Input parameters:
%   pf         Fourier transformed projections are returned from cryo_pft.
%   rotations  The rotation matrix corresponding to each of the
%              projections.
%   max_shift  The maximal 1D shift to look for for each pair of Fourier
%              rays. XXX can be eliminated?
%   shift_step The step size (resolution) which we use to search for 1D
%              shift between Fourier rays. Default: 1 pixel.
%   memoryfactor  If there are N projections, then the system of
%              eqations solved for the shifts is of size 2N x N(N-1)/2 (2N
%              unknowns and N(N-1)/2 equations). This may too big if N is
%              large. If memoryfactor between 0 and 1, then it i the
%              fraction of equation to retain. That is, the system of 
%              equations solved will be of size 2N x N*(N-1)/2*memoryfactor.
%              If memoryfactor is larger than 100, then the number of
%              equations is estimated such the memory used by the equations
%              is roughly memoryfactor megabytes. Default is 1 (use all 
%              equations).
%   shift_2d_ref    Refernce 2D shifts. Used for testing and debugging.
%   verbose       Nonzero to show debug figures and messages. Defualt is 0.
%
% See test_cryo_estimate_shifts for a usage example.
%
% Yoel Shkolnisky, January 2015.

if ~exist('shift_step','var')
    shift_step=1;
end

if ~exist('memoryfactor','var')
    memoryfactor=10000;
end

if ~exist('shifts_2d_ref','var') || isempty(shifts_2d_ref)
    shifts_2d_ref=-1;
end

if ~exist('verbose','var')
    verbose=0;
end

if memoryfactor<0 || (memoryfactor>1 && memoryfactor<100)
    error('subsamplingfactor must be between 0 and 1 or larger than 100');
end

n_theta=size(pf,2);
n_theta2=n_theta/2;
pf=[flipdim(pf(2:end,n_theta/2+1:end,:),1) ; pf(:,1:n_theta/2,:) ];
n_projs=size(rotations,3);

NequationsTotal = ceil(n_projs*(n_projs-1)/2); % Number of equations that will be used to estimation the shifts
memtotal=NequationsTotal*2*n_projs*8; % Estimated memory requirements for the full system of equation.
                                     % This ignores the sparsity of the system, since backslash seems to
                                     % ignore it.

                                     
if memoryfactor<=1                                     
    Nequations = ceil(n_projs*(n_projs-1)*memoryfactor/2); % Number of equations that will be used to estimation the shifts
else
    subsamplingfactor=(memoryfactor*10^6)/memtotal; % By how much we need to 
        % subsample the system of equations in order to use roughly
        % memoryfactor MB.
    if subsamplingfactor<1
        Nequations = ceil(n_projs*(n_projs-1)*subsamplingfactor/2);
    else
        Nequations = NequationsTotal;
    end
end


if Nequations<n_projs
    warning('Too few equations. Increase memoryfactor. Setting Nequations to n_projs');
    Nequations=n_projs;
end

if Nequations<2*n_projs
    warning('Number of equations is small. Consider increase memoryfactor.');
end


% Allocate storage for the equations of determining the 2D shift of each
% projection. The shift equations are represented using a sparse matrix,
% since each row in the system contains four non-zeros (as it involves
% exactly four unknowns). 
% The variables below are used to construct this sparse system. The k'th
% non-zero element of the equations matrix is stored at index 
% (shift_I(k),shift_J(k)).
shift_I=zeros(4*Nequations,1);  % Row index for sparse equations system.
shift_J=zeros(4*Nequations,1);  % Column index for sparse equations system.
shift_eq=zeros(4*Nequations,1); % The coefficients of the center estimation
                                % system ordered as a single vector.                                    
shift_b=zeros(Nequations,1);    % Right hand side of the system.

% Prepare the shift phases
n_shifts=ceil(2*max_shift/shift_step+1); % Number of shifts to try.
rmax=(size(pf,1)-1)/2;
rk=-rmax:rmax; rk=rk(:);
rk2=rk(1:rmax);
shift_phases=zeros(rmax,n_shifts);
for shiftidx=1:n_shifts
    shift=-max_shift+(shiftidx-1)*shift_step;
    shift_phases(:,shiftidx)=exp(-2*pi*sqrt(-1).*rk2.*shift./(2*rmax+1));
end

H=sqrt(abs(rk)).*exp(-rk.^2/(2*(rmax/4).^2)); 
H=H(:);  % Filter for common-line detection.

dtheta=pi/n_theta2; % Not 2*pi/n_theta, since we divided n_theta by 2 
                    % (resulting in n_theta2) to take rays of length
                    % 2*n_r-1. 


% Generate indices of Nequations pairs of common lines.
% Generate all pairs (i,j) such that j>i.
[pairsI,pairsJ]=meshgrid(1:n_projs,1:n_projs);
idxI=pairsI(pairsJ>pairsI);
idxJ=pairsJ(pairsJ>pairsI);
Icl=[idxI,idxJ];
% Pick Nequations indices from I at random.
rp=randperm(size(Icl,1));
Icl=Icl(rp(1:Nequations),:);
%Icl=Icl((1:Nequations),:);

% Iterate over the common lines pairs in Icl and for each pair find the 1D
% relative shift between the two Fourier lines in the pair.
for shift_equation_idx=1:Nequations
    % Find the pair of projections participating in the current common line
    % pair.
    idxi=Icl(shift_equation_idx,1);  % Index of projection i in the pair.
    idxj=Icl(shift_equation_idx,2);  % Index of projection j in the pair.
    
    % Extract the indices of the common line between Pi and Pj.
    Ri=rotations(:,:,idxi);
    Rj=rotations(:,:,idxj);
    [cij,cji]=commonline_R(Ri.',Rj.',n_theta);
    
    % To match cryo_clmatrix, cij is always less than PI and cji may be be
    % larger than PI.
    if cij>=n_theta/2
        cij=cij-n_theta/2;
        cji=cji-n_theta/2;
    end
    if cji<0
        cji=cji+n_theta;
    end
    
    cij=cij+1; cji=cji+1; % Since commonline_R returns zero-based indices.
    
    % Extract the Fourier rays that correspond to the common line.
    pfi=pf(:,cij,idxi);
    
    isPfjFlipped=0;  % Is the common line in image j in the positive
                     % direction of the ray (isPfjflipped=0) or in the
                     % negative direction (isPfjflipped=1).    
    if cji<=n_theta2
        pfj=pf(:,cji,idxj);
    else
        pfj=pf(:,cji-n_theta2,idxj);
        isPfjFlipped=1;
    end
    
    % Find 1D shift between pfi and pfj.
    
    pfi=pfi.*H;
    pfi(rmax:rmax+2,:)=0;
    pfi=cryo_raynormalize(pfi);
    pfi=pfi(1:rmax,:);  % Take half ray plus the DC
        
    pfj=pfj.*H;
    pfj(rmax:rmax+2,:)=0;
    pfj=cryo_raynormalize(pfj);
    pfj=pfj(1:rmax);

    pfi_flipped=conj(pfi);    
    pfi_stack=bsxfun(@times,pfi,shift_phases);
    pfi_flipped_stack=bsxfun(@times,pfi_flipped,shift_phases);
                
    C1=2*real(pfi_stack'*pfj);
    C2=2*real(pfi_flipped_stack'*pfj);
    
    [sval1,sidx1]=max(C1(:));
    [sval2,sidx2]=max(C2(:));
    
    if sval1>sval2 % Rays match in the same orientation.
        dx=-max_shift+(sidx1-1)*shift_step;
    else %Rays match in opposite orientation.
        dx=-max_shift+(sidx2-1)*shift_step;
    end
    
    % Create a shift equation for the projections pair (k1,k2).
    idx=4*(shift_equation_idx-1)+1:4*shift_equation_idx;
    shift_alpha=(cij-1)*dtheta;  % Angle of common ray in projection i.
    shift_beta= (cji-1)*dtheta;  % Angle of common ray in projection j.
    shift_I(idx)=shift_equation_idx; % Row index to construct the sparse equations.
    shift_J(idx)=[2*idxi-1 2*idxi 2*idxj-1 2*idxj]; % Columns of the shift variables that correspond to the current pair (k1,k2).
    shift_b(shift_equation_idx)=dx; % Right hand side of the current equation.
    
    % Compute the coefficients of the current equation.
    if ~isPfjFlipped
        shift_eq(idx)=[sin(shift_alpha) cos(shift_alpha) -sin(shift_beta) -cos(shift_beta)];
    else
        shift_beta=shift_beta-pi; % In the derivation we assume that all 
                        % angles are less than PI where angles larger than
                        % PI are assigned negative orientation.
        shift_eq(idx)=[-sin(shift_alpha) -cos(shift_alpha) -sin(shift_beta) -cos(shift_beta)];
    end
end

shift_equations=sparse(shift_I,shift_J,shift_eq,Nequations,2*n_projs);
shift_equations=[shift_equations shift_b(1:shift_equation_idx)];
%est_shifts=shift_equations(:,1:end-1)\shift_equations(:,end);
est_shifts=lsqr(shift_equations(:,1:end-1),shift_equations(:,end),1.0e-8,100);
est_shifts=full(transpose(reshape(est_shifts,2,n_projs)));

if ~isscalar(shifts_2d_ref) 
    if nnz(shift_equations(:,1:end-1))<10^7
            [~,~,V]=svd(full(shift_equations(:,1:end-1)));
            s1=reshape(shifts_2d_ref.',2*n_projs,1);
            s2=reshape(est_shifts.',2*n_projs,1);
            V=V(:,1:end-3); % Null space of shift_equations.
            % Compute the difference between the true and estimated shifts in
            % the subspace that is orthogonal to the null space of
            % shift_equations.

            if norm(V.'*s1)>1.0e-12
                log_message('cryo_estimate_shifts error: %8.5e',...
                    (norm(V.'*(s1-s2))/norm(V.'*s1)));
            else
                % Print absolute error
                log_message('norm(V.''*s1) = %7.5e',norm(V.'*s2));
            end
    else
        log_message('Not comparing to reference shifts - too many equations\n');
    end
end

if verbose
    if n_projs<=100
        s=svd(full(shift_equations(:,1:end-1)));
        log_message('Singular values of the shift system of equations:');
        log_message('%d  ',fliplr(s.'));
        figure;
        bar(s(end:-1:end-19))
    end
end


