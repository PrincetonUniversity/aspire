function [ clstack,corrstack, shift_equations,shift_equations_map]...
    = commonlines_gaussian_C2(pf,max_shift,shift_step,min_dist_cls)
% Detect the commonlines by searching for the maximun cross-correlation
% between the rays on the polar Fourier Transforms. Gaussian filter applied
% to damp the noise in the high frequency domain.

% Input:
%
%   pf       3D array where each image pf(:,:,k) corresponds to the Fourier
%            transform of projection k.
%   max_shift       Maximal 1D shift (in pixels)  to search between
%       common-lines. Default: 15.
%   shift_step      Resolution of shift estimation in pixels. Note that
%        shift_step can be any positive real number. Default: 1.
%
%
% Output:
%
%   clstack     Common lines matrix. (k1,k2) and (k2,k1) contain the index
%       of the common line of projections k1 and k2. (k1,k2)  contains the
%       index of the common line in the projection k1. (k2,k1) contains the
%       index of the common line in k2.
%   corrstack   The correlation of the common line between projections k1
%       and k2. Since corrstack is symmetric, it contain entries only above
%       the diagonal. corrstack(k1,k2) measures how ''common'' is the between
%       projections k1 and k2. Large value means high-similariry.
%   shift_equations  System of equations for the 2D shifts of the
%       projections. This is a sparse system with 2*n_proj+1 columns. The
%       first 2*n_proj columns correspond to the unknown shifts dx,dy of
%       each projection. The last column is the right-hand-side of the
%       system, which is the relative shift of each pair of common-lines.
%   shift_equations_map   2D array of size n_proj by n_proj. Entry (k1,k2)
%       is the index of the equation (row number) in the array
%       "shift_equations" that corresponds to the common line between
%       projections k1 and k2. shift_map is non-zero only for k1<k2.
%
%
%
% This is a revision on cryo_clmatrix_v3.m by Yoel.
% Added feature: accelerate the searching process by comparing one slice
% with all the other slices simultaneously in each iteration.
%
% Lanhui Wang, July 2, 2013

n_shift=2*max_shift/shift_step+1; % Number of shifts to try.

[n_r,n_theta,n_proj]=size(pf);

% pf is of size n_rxn_theta. Convert pf into an array of size
% (2xn_r-1)xn_theta, that is, take then entire ray through the origin.
% note that the minus one in 2xn_r-1 is because the dc term appers both in
% an angle theta and in angle theta+pi.
% pf columns are (n_theta/2+1,n_theta/2+2,...,n_theta,1,..n_theta/2-1,n_theta/2)
pf=[flipdim(pf(2:end,n_theta/2+1:end,:),1) ; pf(:,1:n_theta/2,:) ];
temp = pf;
% allocate pf to hold a double-cover of all n_theta rays
pf = zeros(2*n_r-1,n_theta,n_proj);
pf(:,1: n_theta/2, :) = temp;
pf(:, n_theta/2+1: end, :) = flipdim(temp, 1); % the second-cover of all rays are simply heading the opposite side
pf = reshape(pf,2*n_r-1,n_theta*n_proj);


%% Apply Gaussian filter
rmax = n_r - 1; %exclude the dc term
rk=-rmax:rmax; rk=rk(:);
H=sqrt(abs(rk)).*exp(-rk.^2/(2*(rmax/4).^2));
pf=bsxfun(@times,pf,H);
% the volume is real, hence fourier-trnsform is conjugate-symmetric. Hence
% all we need are half-rays.
pf=pf(1:rmax,:);
N = rmax;
coefficients=reshape(pf,N,n_theta,n_proj);

%% Allocate variables used for shift estimation

shifts_1d=zeros(n_proj,n_proj);     % Estimated 1D shift between common-lines.

% Based on the estimated common-lines, construct the equations for
% determining the 2D shift of each projection. The shift equations are
% represented using a sparse matrix, since each row in the system contains
% four non-zeros (as it involves exactly four unknowns).
% The variables below are used to construct this sparse system. The k'th
% non-zero element of the equations matrix is stored at index
% (shift_I(k),shift_J(k)).
shift_I=zeros(4*n_proj^2,1);  % Row index for sparse equations system.

shift_J=zeros(4*n_proj^2,1);  % Column index for sparse equations system.

shift_eq=zeros(4*n_proj^2,1); % The coefficients of the center estimation
% system ordered as a single vector.

shift_equations_map=zeros(n_proj); % Entry (k1,k2) is the index of the
% euqation for the common line of projections k1 and k2.

shift_equation_idx=1;  % The equation number we are currently processing.
shift_b=zeros(n_proj*(n_proj-1)/2,1);   % Right hand side of the system.
dtheta=2*pi/n_theta;


%% Search for commonlines
clstack=zeros(n_proj,n_proj,2);      % Common lines-matrix.
corrstack=zeros(n_proj,n_proj,2);    % Correlation coefficient for each common-line.
% normalization
for k=1:n_theta
    for l=1:n_proj
        coefficients(:,k,l)=coefficients(:,k,l)/norm(coefficients(:,k,l));
    end
end
C=coefficients(:,1:n_theta/2,:);
C=reshape(C,rmax,n_theta/2*n_proj); % stack all the Fourier slices
rk2=rk(1:rmax);

% Accelerate the searching process by comparing one slice
% with all the other shifted slices in each iteration.
for k=1:n_proj
    temp_coef=zeros(N,n_theta,n_shift);
    % Generate all the shifted copies of k_th slice.
    for shiftidx=1:n_shift
        shift=-max_shift+(shiftidx-1)*shift_step;
        % a shift in time domain corresponds to rotation in the frequency
        % domain
        shift_phases=exp(-2*pi*sqrt(-1).*rk2.*shift./(2*rmax+1));
        temp_coef(:,:,shiftidx)=bsxfun(@times,coefficients(:,:,k),shift_phases);
    end
    temp_coef=reshape(temp_coef,N,n_theta*n_shift);
    % Compute the cross correlations with all other slices.
    corr=temp_coef'*C(:,n_theta/2*k+1:end);
    % we are comparing slice #k with remaining $n_proj-k$ slices.
    %Every pair of slices induces $n_shift*n_theta^2/2$ correlation scores, placed as a column.
    % the total $$n_shift*n_theta^2$$ can be deduced due to the
    % conjugate-symmetric proprty of the fourier transform (we do so manully)
    corr=reshape(corr,n_shift*n_theta^2/2,n_proj-k);
    
    %GABI : to support C2, we need to pick the second hieghst score as
    %well. but before we want to mask-out close entries
    
    % pick the largest correlation for every pair of slices
    [correlation,idx]=max(real(corr));
    corrstack(k,k+1:n_proj,1)=correlation;
    [cl1,est_shift,cl2]=ind2sub([n_theta n_shift n_theta/2],idx);
    
    maskHole = min_dist_cls*n_theta/360;
    rowMask1 = mod(bsxfun(@plus,cl1',-maskHole:maskHole) -1,n_theta) + 1;
    colMask1 = mod(bsxfun(@plus,cl2',-maskHole:maskHole) -1,n_theta) + 1;
    
    mask     = ones(n_theta,n_shift,n_theta/2,n_proj-k);
    for k_proj=1:(n_proj-k)
            rm1 = rowMask1(k_proj,:);
            cm1 = colMask1(k_proj,:);
            
            cm2 = cm1;
            cm2(cm2<=n_theta/2) = [];
            cm2 = cm2 - n_theta/2;
            
            cm1(cm1>n_theta/2)  = [];
            
            rm2 = mod(rm1+n_theta/2-1,n_theta)+1;
            
            mask(rm1,est_shift(k_proj),cm1,k_proj) = 0;
            mask(rm2,est_shift(k_proj),cm2,k_proj) = 0;
    end
    mask = reshape(mask,n_shift*n_theta^2/2,n_proj-k);
    
    maskedCorr = corr.*mask;
    
    [correlation_g,idx_g]=max(real(maskedCorr));
    corrstack(k,k+1:n_proj,2)=correlation_g;
    [cl1_g,est_shift_g,cl2_g]=ind2sub([n_theta n_shift n_theta/2],idx_g);
    
    % account for the indeces and fourier's conjugate symmetry. Namely:
    % a. cl1 reported indeces actually correspond to n_theta values :
    %    (n_theta/2+1,n_theta/2+1,...,n_theta,1,..n_theta/2-1,n_theta/2)
    % b. cl2 reported indeces actually correspond to n_theta/2 values :
    %    (n_theta/2+1,n_theta/2+1,...,n_theta)
    %     cl2(cl1>n_theta/2)=cl2(cl1>n_theta/2)+n_theta/2;
    %     cl1(cl1>n_theta/2)=cl1(cl1>n_theta/2)-n_theta/2;
    %
    %     cl2_g(cl1_g>n_theta/2)=cl2_g(cl1_g>n_theta/2)+n_theta/2;
    %     cl1_g(cl1_g>n_theta/2)=cl1_g(cl1_g>n_theta/2)-n_theta/2;
    
    shifts_1d(k,k+1:n_proj,1)=-max_shift+(est_shift-1)*shift_step;
    shifts_1d(k,k+1:n_proj,2)=-max_shift+(est_shift_g-1)*shift_step;
    
    clstack(k,k+1:n_proj,1)=cl1;
    clstack(k+1:n_proj,k,1)=cl2';
    
    clstack(k,k+1:n_proj,2)=cl1_g;
    clstack(k+1:n_proj,k,2)=cl2_g';
end

for k1=1:n_proj
    for k2=k1+1:n_proj
        idx=4*(shift_equation_idx-1)+1:4*shift_equation_idx;
        cl1=clstack(k1,k2);
        cl2=clstack(k2,k1);
        shift_alpha=(cl1-1)*dtheta;  % Angle of common ray in projection 1.
        shift_beta= (cl2-1)*dtheta;  % Angle of common ray in projection 2.
        shift_I(idx)=shift_equation_idx; % Row index to construct the sparse equations.
        shift_J(idx)=[2*k1-1 2*k1 2*k2-1 2*k2]; % Columns of the shift variables that correspond to the current pair (k1,k2).
        shift_b(shift_equation_idx)=shifts_1d(k1,k2); % Right hand side of the current equation.
        
        % Compute the coefficients of the current equation.
        if shift_beta<pi-1e-13
            shift_eq(idx)=[sin(shift_alpha) cos(shift_alpha) -sin(shift_beta) -cos(shift_beta)];
        else
            shift_beta=shift_beta-pi; % In the derivation we assume that all angles are less
            % than PI where angles larger than PI are assigned
            % nigative orientation.
            shift_eq(idx)=[-sin(shift_alpha) -cos(shift_alpha) -sin(shift_beta) -cos(shift_beta)];
        end
        
        shift_equations_map(k1,k2)=shift_equation_idx;  % For each pair (k1,k2), store the index of its equation.
        shift_equation_idx=shift_equation_idx+1;
    end
end

shift_equation_idx=shift_equation_idx-1;
shift_equations=sparse(shift_I(1:4*shift_equation_idx),...
    shift_J(1:4*shift_equation_idx),shift_eq(1:4*shift_equation_idx),...
    shift_equation_idx,2*n_proj);
shift_equations=[shift_equations shift_b(1:shift_equation_idx)];