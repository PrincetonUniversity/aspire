function [est_shifts,shift_equations]=...
            cryo_estimate_shifts_ref(pf,rotations,clstack,...
                                            max_shift,shift_step)
%
% CRYO_ESTIMATE_SHIFTS_REF  Esimate 2D shifts in projections.
%
% [est_shifts,shift_equations]=...
%            cryo_estimate_shifts_ref(pf,rotations,clstack,...
%                                            max_shift,shift_step)
%   Estimate 2D shifts in projections using the estimated rotations and the
%   common line matrix.
%   The function processes the projections exactly as cryo_clmatrix_gpu and
%   returns shift equations which are identical to those returned by
%   cryo_clmatrix_gpu. This code is used as a reference code used for the
%   development of cryo_estimate_shifts.


n_theta=size(pf,2);
n_theta2=n_theta/2;
pf=[flipdim(pf(2:end,n_theta/2+1:end,:),1) ; pf(:,1:n_theta/2,:) ];
n_projs=size(rotations,3);
%Nequations = min(n_projs*(n_projs-1)/2,10*n_projs); % Number of equations that will be used to estimation the shifts
Nequations = n_projs*(n_projs-1)/2;

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

% Find "more reliable" common lines. 
% We find the 10*K pairs of common lines that are most reliable. These will
% be used to form a K x 10K system of equations for the 2D shift each each
% projection. We use only 10K pair of common lines and not all K choose 2
% pairs so we can resolve for shifts even when K is large (for otherwise
% the system of equaitons is too large to be solved).

clerr=syncconsistency(rotations,clstack,n_theta);
%[sortedclerr,I]=sort(clerr);
sortedclerr=clerr;
I=1:numel(clerr);
Icl=I(1:Nequations); % Indices of common lines pairs for which we will for equations.

% Iterate over the common lines pairs in Icl and for each pair find the 1D
% relative shift between the two Fourier lines in the pair.
for shift_equation_idx=1:numel(Icl)
    % Find the pair of projections participating in the current common line
    % pair.
    idxi=sortedclerr(shift_equation_idx,1);  % Index of projection i in the pair.
    idxj=sortedclerr(shift_equation_idx,2);  % Index of projection j in the pair.
    
    % Extract the indices of the common line between Pi and Pj.
    cij=clstack(idxi,idxj); % Common line between Pi and Pj in Pi.
    cji=clstack(idxj,idxi); % Common line between Pi and Pj in Pj.
    
    % Extract the Fourier rays that correspond to the common line.
    isPfiFlipped=0; % Is the common line in image i in the positive 
                    % direction of the ray (isPfiflipped=0) or in the
                    % negative direction (isPfiflipped=1).
    if cij<=n_theta2
        pfi=pf(:,cij,idxi);
    else
        pfi=pf(:,cij-n_theta2,idxi);
        %pfi=conj(pfi);
        isPfiFlipped=1;
    end

    isPfjFlipped=0;
    if cji<=n_theta2
        pfj=pf(:,cji,idxj);
    else
        pfj=pf(:,cji-n_theta2,idxj);
        %pfj=conj(pfj);
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
    if ~xor(isPfiFlipped,isPfjFlipped)
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
est_shifts=shift_equations(:,1:end-1)\shift_equations(:,end);
est_shifts=transpose(reshape(est_shifts,2,n_projs));

if n_projs<=100
    s=svd(full(shift_equations(:,1:end-1)));
    log_message('Singular values of the shift system of equations:');
    log_message('%d  ',fliplr(s.'));
    figure;
    bar(s(end:-1:end-19))
end
