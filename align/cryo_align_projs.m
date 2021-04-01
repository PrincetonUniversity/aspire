function [Rots_est,Shifts,corrs,err_Rots,err_Shifts] = cryo_align_projs(sym,n_sym,projs,vol,N_ref,isshift,G,true_Rots,true_Rots_J,true_Shifts,Rots,verbose)
%% This function aligns given projection in a given volume.
%% input: 
%{
sym- symmetry tipe- 'C'/ 'D'/ 'T'/ 'O'/ 'I'.
n_s- symmetry value. for cubic symmetries it doesn't metter (don't have 
     to enter).
projs- projection images to align in a given volume.
vol- reference volume for the alignment.
N_ref- number of reference projections for the alignment. (default is 50).
isshift- enter 1 if you wish to estimate the translations of the projections. 
         and enter 0 if not. 


G- size=(3,3,n) all n symmetry group elemnts. 
true_Rots- the true rotations of projs.
true_Rots_J- the true rotations of projs in case of reflection.
true_Shifts- the true shifts-(dx,dy) of projs.
Rots - size (3x3xsz_Rots), is a set of candidate rotations. (don't have to
       enter).
verbose- enter some number different from 0 if you wish to print log 
       messages. (default is 0).
%}
%% output:
%{
Rots_est- the rotation matrices of the alingment of projs (3x3xsize(projs,3)).
Shifts- (size(projs,3)x2) the 2D shift of each projection, first column
       contained the shift in the x-axis, and the secound colunm in the y-axis.
corrs- (size(projs,3)x2) the i'th entry of the first column contains 
       the correlation of the common lines between the i'th image and all the 
       reference images induced by the best matching rotation. The i'th entry of 
       the second column contains the mean matching correlation over all tested 
       rotations.
candidate_rots- set of candidate rotations that from them the estimated
       rotation for each projection is chosen. 
err_Rots- error calculation between the true rotations and the estimated
       rotations.
err_Shifts- error calculation between the true shifts and the estimated
       shifts, in x and y axis.
%}
%% Check parameters:
if ~exist('N_ref','var') || isempty(N_ref)
    N_ref = 30;
end

er_calc = 1;
if sym == 'C' && n_sym == 1
    er_calc = 1;
    G = eye(3);
elseif ~exist('G','var') || isempty(G)
    er_calc = 0;
end
ref_true_rot = 1;
if ~exist('true_Rots','var') || isempty(true_Rots)
    ref_true_rot = 0;
end      
ref_true_rot_J = 1;
if ~exist('true_Rots_J','var') || isempty(true_Rots_J)
    ref_true_rot_J = 0;
end 
refshift = 1;
if ~exist('true_Shifts','var') || isempty(true_Shifts)
    refshift = 0;
end
can_Rots = 1;
if ~exist('Rots','var') || isempty(Rots)
    can_Rots = 0;
end
if ~exist('verbose','var') || isempty(verbose)
    verbose = 0;
end
currentsilentmode = log_silent(verbose == 0);

n = size(vol,1);
L = 360;

%% Compute polar Fourier transform of projs:
n_r = ceil(n/2);
log_message('Start computing polar Fourier transforms of input projections Using n_r=%d L=%d.',n_r,L);
projs_hat = cryo_pft(projs,n_r,L,'single');
log_message('Computing polar Fourier transform done');

% Normalize polar Fourier transforms
log_message('Start normalizing polar Fourier transform of input projections (cryo_raynormalize)');
projs_hat = cryo_raynormalize(projs_hat);
log_message('Normalizing done');
n_projs = size(projs_hat,3);

%% Generate candidate rotations and reference projections:
log_message('Generating %d reference projections.',N_ref);

if can_Rots == 0
    Rots = genRotationsGrid(75);
end

candidate_rots = Rots;
N_rot = size(candidate_rots,3);
log_message('Using %d candidate rotations for alignment.',N_rot);

rots_ref = Rots(:,:,randperm(N_rot,N_ref));   

ref_projs = cryo_project(vol,rots_ref,n);
ref_projs = permute(ref_projs,[2 1 3]);
rots_ref = permute(rots_ref,[2,1,3]);         % the true rotations.

%% Compute polar Fourier transform of reference projections:
log_message('Start computing polar Fourier transforms of reference projections. Using n_r=%d L=%d.',n_r,L);
refprojs_hat = cryo_pft(ref_projs,n_r,L,'single');
log_message('Computing polar Fourier transform done');

% Normalize polar Fourier transforms
log_message('Start normalizing polar Fourier transform of reference projections (cryo_raynormalize)');
refprojs_hat = cryo_raynormalize(refprojs_hat);
log_message('Normalizing done');


%% Compute the common lines between the candidate rotations and the reference rotations:
Ckj = (-1)*ones(N_rot,N_ref);   % In the coordinates of candidate_rots.
Cjk = (-1)*ones(N_rot,N_ref);   % In the coordinates of rots_ref.
Mkj = zeros(N_rot,N_ref);       % Pairs of rotations that are not "too close"
for k = 1:N_rot 
    Rk = candidate_rots(:,:,k).';
    for j = 1:N_ref        
        Rj = rots_ref(:,:,j).';
         if sum(Rk(:,3).*Rj(:,3)) < 0.999
             [ckj,cjk] = commonline_R2(Rk,Rj,L); 
             %%% Convert the returned indices ckj and cjk into 1-based.
             ckj = ckj+1; cjk = cjk+1;         
             Ckj(k,j) = ckj;
             Cjk(k,j) = cjk;
             Mkj(k,j) = 1;
         end    
    end
end

%% Generate shift grid:
% generating a shift grid on the common lines, and choosing the shift
% that brings the best correlation in the comparisson between the common 
% lines. 
% after applying polar FFT on each projection, the shift in quartesian
% coordinates- (delta(x),delta(y)) becomes a shift only in the r variable
% in the common lines (the common line have a specific theta so we have to
% consider a shift only in the r variable.
% the equation for the shift phase in the common lines is:
% exp((-2*pi*i)*r*delta(r)).
max_s = round((0.2)*size(projs_hat,1));     % set the maximum shift.
s_step = 0.5;
n_shifts = (2/s_step)*max_s + 1;            % always odd number (to have zero value without shift).
max_r = size(projs_hat,1);
%log_message('Using max_shift=%d  shift_step=%5.1e, n_shifts=%d for projections of size_r=%d',max_s,s_step,n_shifts,max_r);
s_vec = linspace(-max_s,max_s,n_shifts);    % those are the shifts in the r variable in the common lines.  
r_vec = (0:max_r-1);
s_phases = exp(-2*pi*sqrt(-1).*r_vec'*s_vec./(2*max_r+1));    % size of (n_rXn_shift)

%% Main loop- compute the cross correlation: 
% computing the correlation between the common line, first choose the best
% shift, and then chose the best rotation.
Rots_est = zeros(3,3,n_projs);
corrs = zeros(n_projs,2);                   % Statistics on common-lines matching.
Shifts = zeros(2,n_projs);
dtheta = 2*pi/L;
if ref_true_rot ~= 0 || ref_true_rot_J ~= 0
    err_Rots = zeros(n_projs,1);
end
if refshift ~= 0 
    err_Shifts = zeros(2,n_projs);
end
for projidx = 1:n_projs
    cross_corr_m = zeros(N_rot,N_ref);
    for j = 1:N_ref
        iidx = find(Mkj(:,j)~=0);
        conv_hat = bsxfun(@times,conj(projs_hat(:,Ckj(iidx,j),projidx)),refprojs_hat(:,Cjk(iidx,j),j));   % size of (n_rxsize(iidx))
        temp_corr = real(s_phases'*conv_hat);      % size of (n_shiftXsize(iidx)).
        cross_corr_m(iidx,j) = max(temp_corr).';
    end
    %%% calculating the mean of each row in cross_corr_m:
    cross_corr = sum(cross_corr_m,2)./sum(cross_corr_m>0,2); 
    
    %% Find estimated rotation:
    [bestRscore,bestRidx] = max(cross_corr);
    meanRscore = mean(cross_corr);
    corrs(projidx,1) = bestRscore;
    corrs(projidx,2) = meanRscore;
    Rots_est(:,:,projidx) = candidate_rots(:,:,bestRidx);

    %%% Error calc for estimated rotation:
    if ref_true_rot ~= 0 && er_calc ~= 0
        g_est_t = Rots_est(:,:,projidx)*true_Rots(:,:,projidx).';
        n_g = size(G,3);
        dist = zeros(n_g,1);
        for g_idx = 1:n_g
            dist(g_idx,1) = norm(g_est_t-G(:,:,g_idx),'fro');
        end
        [~,min_idx] = min(dist);
        g_est = G(:,:,min_idx);

        R_est = g_est.'*Rots_est(:,:,projidx);
        R = true_Rots(:,:,projidx)*(R_est.');
        err_Rots(projidx,:) = (acosd((trace(R)-1)/2));     %in deg.  
    end
        %%% Error calculation for reflection case:
        % if there is a reflection between the projection and the volume
        % then, the relation is R_est=gJRJ.
    if ref_true_rot_J ~= 0 && er_calc ~= 0
        J3 = diag([1 1 -1]);
        g_est_t = Rots_est(:,:,projidx)*(J3*true_Rots_J(:,:,projidx)*J3).';
        n_g = size(G,3);
        dist = zeros(n_g,1);
        for g_idx = 1:n_g
            dist(g_idx,1) = norm(g_est_t-G(:,:,g_idx),'fro');
        end
        [~,min_idx] = min(dist);
        g_est = G(:,:,min_idx);

        R_est = J3*g_est.'*Rots_est(:,:,projidx)*J3;
        R = true_Rots_J(:,:,projidx)*(R_est.');
        err_2 = (acosd((trace(R)-1)/2));     %in deg.  
        if err_Rots(projidx,:) ~= 0
            err_Rots(projidx,:) = min(err_Rots(projidx,:),err_2);
        else
            err_Rots(projidx,:) = err_2; 
        end
    end      
               
    

    %% Find estimated shift:
    % by least-squares on the estimated rotation with the reference projections. 
    if isshift ~= 0 
        idx = find(Mkj(bestRidx,:)==1);
        [n,~] = size(idx);
        shift_eq = zeros(n,2);
        shift = zeros(n,1);
        i = 1;
        for j = idx
            conv_hat = bsxfun(@times,conj(projs_hat(:,Ckj(bestRidx,j),projidx)),refprojs_hat(:,Cjk(bestRidx,j),j));
            temp_corr = real(s_phases'*conv_hat);      % size of (n_shiftX1).
            [~,I] = max(temp_corr); 
            s_idx = I;
            shift(i,1) = s_vec(1,s_idx);
            theta = (Ckj(bestRidx,j)-1)*dtheta;     % Angle of the common line.
            shift_eq(i,1) = sin(theta);
            shift_eq(i,2) = cos(theta);
            i = i+1;
        end

        Shifts(:,projidx) = shift_eq\shift;
        %%% Error calc for estimated shifts:
        if refshift ~= 0 
            err_Shifts(1,projidx) = norm(true_Shifts(projidx,1)-Shifts(1,projidx),2);
            err_Shifts(2,projidx) = norm(true_Shifts(projidx,2)-Shifts(2,projidx),2);
        end
    end
end
if ref_true_rot ~= 0 && er_calc ~= 0
    mean_err = mean(err_Rots);
    log_message('Mean error in estimating the rotations of the projections is %5.3f in degrees.',mean_err);
end
if isshift ~= 0 && refshift ~= 0
    mean_err_shift = mean(mean(err_Shifts));
    log_message('Mean error in estimating the translations of the projections is %5.3f.',mean_err_shift);
end

log_silent(currentsilentmode);
end


           



    
