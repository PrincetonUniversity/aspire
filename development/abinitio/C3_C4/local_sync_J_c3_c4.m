function [vijs,viis,im_inds_to_remove,pairwise_inds_to_remove,npf,projs,refq,ref_shifts] = local_sync_J_c3_c4(n_symm,Rijs,Riis,npf,...
                                projs,is_remove_non_rank1,remov_percent,refq,ref_shifts)

% Local J-synchronization of all relative orientations.
%
% Input parameters:
%   n_symm               Either 3 (for c_3) or 4 (for c_4)
%   Rijs                 A 3x3xn_choose_2 array holding the estimates for the
%                        relative orientations. Specifically, each slice Rij
%                        is equal to either $Ri^{T}g^{sij}Rj$, or $JRi^{T}g^{sij}RjJ$
%                        where for n_symm=3 sij is either 0,1 or 2 and for n_symm=4 sij is either 0,1,2 or 3 
%
%   Riis                 A 3x3xn array holding the estimates for the
%                        self relative orientations. Specifically, each slice Rii
%                        is equal to either $Ri^{T}g^{si}Rj$, or $J Ri^{T}g^{si}Rj J$
%                        where for n_symm=3 each si equals 1 or 2, and for
%                        n_symm=4 each si equals 1 or 3
%   npf                  (Optional, required if is_remove_non_rank1==1) 
%                        A 3D array where each image npf(:,:,i) corresponds to the Fourier
%                        transform of projection i.
%   projs                (Optional, required if is_remove_non_rank1==1) 
%                        A 3D array of projection images
%   is_remove_non_rank1  (Optional) zero or one. Whether to remove or not images that
%                        induce too many non-rank-1 matrices. Defualt=1
%   remov_percent
%   refq                 (Optional) A 4-by-n table. The i-th column represent the quaternion of
%                        that corresponds to the rotation matrix of the i-th image
%
%   ref_shifts           (Optional) An nx2 table of shifts that each image
%                        has undergone
%
% Output parameters:
%   vijs           A 3x3xn_choose_2 array where each slice holds an estimate for
%                  the corresponding outer-product vi*vj^{T} between the
%                  third rows of matrices Ri and Rj. Each such estimate
%                  might have a spurious J independently of other estimates
%   viis           A 3x3xn array where the i-th slice holds an estimate for
%                  the outer-product vi*vi^{T} between the
%                  third row of matrix Ri with itself. Each such estimate
%                  might have a spurious J independently of other estimates
%  im_inds_to_remove The image indexes that are removed since they induce
%                    too many non rank-1 matrices
%  pairwise_inds_to_remove 
%   npf            Only if provided in input. Returns the input npf where
%                  all images that correspond to image indexes to remove are removed
%   projs          Only if provided in input. Returns the input projs where
%                  all projections that correspond to image indexes to remove are removed
%   refq           Only if provided in input. Returns the input refq where
%                  all queternions that correspond to image indexes to remove are removed
%   ref_shifts     Only if provided in input. Returns the input ref_shifts where
%                  all shifts that correspond to image indexes to remove are removed

log_message('Local J synchronization');

if n_symm ~= 3 && n_symm ~= 4
    error('n_symm may be either 3 or 4');
end

if exist('is_remove_non_rank1','var') && is_remove_non_rank1 == 1
    if ~exist('npf','var')
        error('variable npf must be given if is_remove_non_rank1==true');
    end
    
    if ~exist('projs','var')
        error('variable projs must be given if is_remove_non_rank1==true');
    end
end

if ~exist('is_remove_non_rank1','var')
    is_remove_non_rank1 = true;
end

if ~exist('non_rank1_remov_percent','var')
    remov_percent = 0.25;
end

nImages = size(Riis,3);
assert(nchoosek(nImages,2) == size(Rijs,3));

nrank1 = 0;
e1 = [1 0 0].';
msg = [];

viis = zeros(3,3,nImages);
if n_symm == 3
    % no matter whether Rii=RigRi or Rii=Rig^{2}Ri (and, possibly
    % also J-conjugated), the sum Rii + Rii.' + eye(3) is rank-1 and is equal to the
    % outor product vi*vj^{T} (or J*vi*vj^{T}*J) of third row of matrices Ri
    % and Rj
    for i=1:nImages
        Rii = Riis(:,:,i);
        viis(:,:,i) = (Rii + Rii.' + eye(3))/3;
    end
else % i.e., n_symm = 4
    % no matter whether Rii=RigRi or Rii=Rig^{3}Ri (and, possibly
    % also J-conjugated), the sum Rii + Rii.' is rank-1 and is equal to the
    % outor product vi*vj^{T} (or J*vi*vj^{T}*J) of third row of matrices Ri
    % and Rj
    for i=1:nImages
        Rii = Riis(:,:,i);
        viis(:,:,i) = 0.5*(Rii + Rii.');
    end
end

vijs         = zeros(3,3,nchoosek(nImages,2));
isRank1_ijs  = zeros(nImages,nImages);
stats        = zeros(1,nchoosek(nImages,2));
J = diag([1 1 -1]); % Reflection matrix
TOL = 1.0E-1; % tollerance value
for i=1:nImages
    
    for j=i+1:nImages
        
        t1 = clock;
        ind = uppertri_ijtoind(i,j,nImages);
        
        Rij = Rijs(:,:,ind);
        
        Rii = Riis(:,:,i);
        Rjj = Riis(:,:,j);
        
        % A rank-1 matrix is attained if and only if the
        % following two conditions are satisfied:
        % 1. matrices Rij,Rii,Rjj either are all J-conjugated or none are at all
        % 2. either (Rii = RigRi and Rjj = RjgRj)
        %        or (Rii = Rig^{3}Ri and Rjj = Rjg^{3}Rj).
        
        % There are 2*2*2 possibilities to check, only one of which is of rank one
        % 1. J-conjugate Rii or not
        % 2. J-conjugate Rjj or not
        % 3. transpose Rjj (note that (RjgRj)^T = Rjg^{3}Rj).
        JRiiJ = J*Rii*J;
        c_Rii = {Rii, Rii, Rii, Rii, JRiiJ, JRiiJ, JRiiJ, JRiiJ};
        
        JRjjJ       = J*Rjj*J;
        Rjj_T       = Rjj.';
        JRjjJ_T = J*Rjj_T*J;
        
        c_Rjj = {Rjj, JRjjJ, Rjj_T, JRjjJ_T, Rjj, JRjjJ, Rjj_T, JRjjJ_T};
        
        vij_cands   = zeros(3,3,8);
        svlaues     = zeros(3,8); % three singular values for each possibility
        is_rank1    = false(1,8);
        score_rank1 = zeros(1,8);
        for s = 1:8
            Rii_cand = c_Rii{s};
            Rjj_cand = c_Rjj{s};
            
            if n_symm == 3
                vij_cand = (Rij + Rii_cand*Rij*Rjj_cand + Rii_cand^2*Rij*Rjj_cand^2)/3;
            else
                vij_cand = (Rij + Rii_cand*Rij*Rjj_cand)/2;
            end
            
            vij_cands(:,:,s) = vij_cand;
            svals = svd(vij_cand);
            % meassure how close are the singular values to (1,0,0)
            is_rank1(s) = abs(svals(1)-1) < TOL && sum(abs(svals(2:3)))/2 < TOL;
            score_rank1(s) = norm(svals-e1,2);
            svlaues(:,s) = svals;
        end
        
        if any(is_rank1 == true)
            isRank1_ijs(i,j) = 1;
            nrank1 = nrank1 + 1; % just for stats puproses
        end
        % even if none is rank-1 we still need to choose the best one
        [~,ii] = min(score_rank1);
        vijs(:,:,ind) = vij_cands(:,:,ii);
        stats(ind) = ii;
        
        %%%%%%%%%%%%%%%%%%% debug code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        t2 = clock;
        t = etime(t2,t1);
        bs = char(repmat(8,1,numel(msg)));
        fprintf('%s',bs);
        msg = sprintf('k1=%3d/%3d k2=%3d/%3d t=%7.5f',i,nImages,j,nImages,t);
        fprintf('%s',msg);
        %%%%%%%%%%%%%%%%%%% end of debug code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end
fprintf('\n');

stats_dist = histc(stats,1:8)/numel(stats);
log_message('percentage of rank-1 matrices= %.2f%%', nrank1/nchoosek(nImages,2)*100);
log_message('inner_sync_dist=[%.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f]', stats_dist);


if is_remove_non_rank1
    
    isRank1_ijs = isRank1_ijs + isRank1_ijs.';
    
    % find the image indeces whose relative orientations that they induce
    % involve the largest number of matrices that are not rank 1
    [~,inds_err_rank1] = sort(sum(isRank1_ijs),'ascend');
    
    nRemove = floor(remov_percent*nImages);
    im_inds_to_remove = inds_err_rank1(1:nRemove);
    n = numel(im_inds_to_remove);
    log_message('Removing %d out of %d images that induce non rank-1 matrices',...
        numel(im_inds_to_remove),nImages);
    
    % precalculate the number of relative orientations to remove (inclusion exclusion).
    %Each image to remove is part of nImages-1 relative orientations, but we need to
    % subtract the number of relative orientations that comprise of two images
    % that need to be removed.
    sz = n*(nImages-1)-nchoosek(n,2);
    
    pairwise_inds_to_remove = zeros(1,sz);
    k=1;
    for i=im_inds_to_remove        
        for j=1:nImages
            if j<i
                ind = uppertri_ijtoind(j,i,nImages);
            elseif j==i
                continue;
            else
                ind = uppertri_ijtoind(i,j,nImages);
            end
            pairwise_inds_to_remove(k) = ind;
            k = k+1;
        end
    end
    
    npf(:,:,im_inds_to_remove) = [];
    projs(:,:,im_inds_to_remove) = [];
    viis(:,:,im_inds_to_remove) = [];
    vijs(:,:,pairwise_inds_to_remove) = [];
    
    if exist('refq','var') && ~isempty(refq)
        assert(size(refq,2) == nImages);
        refq(:,im_inds_to_remove) = [];
    end
    
    if exist('ref_shifts','var') && ~isempty(ref_shifts)
        assert(size(ref_shifts,1) == nImages);
        ref_shifts(im_inds_to_remove,:) = [];
    end
else
    % nothing to remove
    im_inds_to_remove = [];
    pairwise_inds_to_remove = [];
end

% if ~params.real_data && params.debug && isfield(params,'inds_err_Rii')
%
%     inds_err_Rii = params.inds_err_Rii;
%     top_inds_err_Rii = inds_err_Rii(1:floor(0.2*numel(inds_err_Rii)));
%
%     [~,locs] = ismember(top_inds_err_Rii,im_inds_to_remove);
%     figure; scatter(locs,zeros(1,numel(locs)));
%
%     %     [~,locs] = ismember(top_inds_err_Rii,inds_err_rank1);
%     %     figure; scatter(locs,zeros(1,numel(locs)));
%
% end
end