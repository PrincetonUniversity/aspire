function [Rijs_out,Rijgs_out] = local_sync_J(Rijs,Rijgs,nImages)

% Local J-synchronization of all relative orientations.
%
% Input parameters:
%   Rijs                 A 3x3xn_choose_2 array holding the estimates for the
%                        relative orientations. Specifically, each slice Rij
%                        is equal to either $Ri^{T}g^{sij}Rj$, or $JRi^{T}g^{sij}RjJ$
%                        where and sij is either 0,1,2 or 3
%
%   Riis                 A 3x3xn array holding the estimates for the
%                        self relative orientations. Specifically, each slice Rii
%                        is equal to either $Ri^{T}g^{si}Rj$, or $J Ri^{T}g^{si}Rj J$
%                        where each si equals 1 or 3.
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

assert(nchoosek(nImages,2) == size(Rijs,3));
assert(nchoosek(nImages,2) == size(Rijgs,3));

Rijs_out  = zeros(size(Rijs));
Rijgs_out = zeros(size(Rijgs));


nrank1 = 0;
e1 = [1 0 0].';
msg = [];

isRank1_ijs  = zeros(nImages,nImages);
stats        = zeros(1,nchoosek(nImages,2));
J = diag([1 1 -1]); % Reflection matrix
TOL = 1.0E-1; % tollerance value
for i=1:nImages
    
    for j=i+1:nImages
        
        t1 = clock;
        ind = uppertri_ijtoind(i,j,nImages);
        
        Rij  = Rijs(:,:,ind);
        Rijg = Rijgs(:,:,ind);
        % Rij+Rij_g must be rank 1. If not J-conjuage either of them
        vij  = (Rij     + Rijg)/2; %should be rank 1
        vijJ = (J*Rij*J + Rijg)/2;
        
       
        % we are only interested in the singular values
        s  = svd(vij);
        sJ = svd(vijJ);
        if (abs(s(1)-1)> TOL || sum(abs(s(2:3)))/2 > TOL) && ...
                (abs(sJ(1)-1)> TOL || sum( abs(sJ(2:3)) )/2 > TOL)
            isRank1_ijs(i,j) = 0;
        else
            isRank1_ijs(i,j) = 1;
            nrank1 = nrank1 + 1; % just for stats puproses
        end
        
        if norm(s-e1,2) < norm(sJ-e1,2)
            Rijs_out(:,:,ind)  = Rij;
            Rijgs_out(:,:,ind) = Rijg;
            stats(ind)   = 0;
        else
            Rijs_out(:,:,ind)  = Rij;
            Rijgs_out(:,:,ind) = J*Rijg*J;
            stats(ind)   = 1;
        end
               
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

stats_dist = histc(stats,0:1)/numel(stats);
log_message('percentage of rank-1 matrices= %.2f%%', nrank1/nchoosek(nImages,2)*100);
log_message('inner_sync_dist=[%.2f %.2f]',stats_dist);

end