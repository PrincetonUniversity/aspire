function cache_filename = cryo_TO_create_cache(symmetry)

% Creates a mat file containing:
%       R - candidate rotations set. 
%       n_theta
%       l_ij_ji_ind - set of common lines in linear indices.
%       l_self_ind - set of self common lines in linear indices.
%
%   Input:
%       symmetry -      either 'T' or 'O'.
%
%   Output:
%       cache_filename - a mat file name containing all required data.
%       
%   Written by Adi Shasha January 2021. 

[n_theta, ~, ~, resolution, viewing_angle, inplane_rot_degree] = cryo_TO_configuration();
[gR, n_gR, n_scl_pairs] = cryo_TO_group_elements(symmetry);

log_message('Generating candidate rotations set');
R = TO_candidate_rotations_set(gR, resolution, viewing_angle, inplane_rot_degree);
n_candidates = size(R,3);

log_message('Computing common lines and self common lines indices sets');
[l_self_ind, l_ij_ji_ind] = TO_cl_scl_inds(gR, n_gR, n_scl_pairs, n_theta, R);

cache_filename = sprintf('%c_symmetry_%d_candidates_cache_test.mat',symmetry, n_candidates);
log_message('Cache file name: %s', cache_filename);
save(cache_filename,'l_self_ind','l_ij_ji_ind','-v7.3','-nocompression');
save(cache_filename,'R','n_theta','-append');    


function R = TO_candidate_rotations_set(gR, resolution, viewing_angle, inplane_rot_degree)

% generates approximatly equally spaced rotations, and then filters rotation 
% with close viewing angle and close inplane rotation according to symmetry
% group elements.
%
%   Input:
%       gR - symmetry group elements: 3 x 3 x group order.
%       resolution - the number of samples per 2*pi.
%                    see genRotationsGrid for more details.
%       viewing_angle - the viewing angle threshold
%       inplane_rot_degree - the inplane rotation degree threshold
%
%   Output:
%       R - set of candidate roation matrices (3 x 3 x n_R) after filtering.
%       
%   Written by Adi Shasha January 2021.

n_gR = size(gR,3);

candidates_set = genRotationsGrid(resolution);
n_candidates = size(candidates_set,3);

close_idx = zeros(1,n_candidates);

for r_i = 1:(n_candidates-1)
    if close_idx(1,r_i) == 1
        continue
    end
    Ri = candidates_set(:,:,r_i);
    for r_j = (r_i+1):n_candidates
        if close_idx(1,r_j) == 1
            continue
        end
        Rj = candidates_set(:,:,r_j);      
        for k = 1:n_gR
            gRj = gR(:,:,k)*Rj;
            if sum(Ri(:,3).*gRj(:,3)) > viewing_angle               % the viewing angles.
                R_inplane = Ri.'*gRj;
                theta = abs(atand(R_inplane(2,1)/R_inplane(1,1)));
                if theta < inplane_rot_degree                       % the inplane rotation.
                    close_idx(1,r_j) = 1;
                end
            end
        end
    end
end

left_idx = close_idx==0;
R = candidates_set(:,:,left_idx);


function [l_self_ind, l_ij_ji_ind] = TO_cl_scl_inds(gR, n_gR, n_scl_pairs, n_theta, R)

% Computes the set of common lines induced by all rotation pairs from R, 
% and the set of self common lines induces by each rotation from R.
%
% Note: computing the set of self common lines requires a certain order of 
% the symmetry group elements.
%
%   Input:
%       gR -            symmetry group elements: 3 x 3 x group order.    
%       n_gR -          group order, also the number of common lines.
%       n_scl_pairs -   number of self common lines.
%       n_theta -       how many radial lines.
%       R -             a set of roation matrices: 3 x 3 x n_R .
%
%   Output:
%       l_self_ind -    set of self common lines in linear indices: 
%                       n_candidates x n_scl_pairs
%       l_ij_ji_ind -   set of common lines in linear indices:
%                       n_candidates_pairs x n_gR
%       
%   Written by Adi Shasha January 2021. 

n_R = size(R, 3);

% array of linear indices of common lines and self common lines
l_self  = zeros(2,n_scl_pairs,n_R);
l_ij_ji = zeros(2,n_gR,n_R*n_R);

for i=1:n_R 
    % Compute self common lines of Ri.
    for k=1:n_scl_pairs
        [l_self(1,k,i),l_self(2,k,i)] = commonline_R(R(:,:,i),R(:,:,i)*gR(:,:,2*k),n_theta);
    end
    for j=1:n_R 
        if i == j
            continue
        end        
        % Compute the common lines induced by rotation matrices Ri and Rj.
        ind = sub2ind([n_R,n_R],i,j);
        for k=1:n_gR
            [l_ij_ji(1,k,ind),l_ij_ji(2,k,ind)] = commonline_R(R(:,:,i),R(:,:,j)*gR(:,:,k),n_theta);
        end
    end
end

l_self  = l_self + 1;
l_ij_ji = l_ij_ji + 1;

% convert to linear indices
l_self_ind  = sub2ind([n_theta,n_theta],l_self(1,:,:),l_self(2,:,:));
l_ij_ji_ind = sub2ind([n_theta,n_theta],l_ij_ji(1,:,:),l_ij_ji(2,:,:));

% reshape
l_self_ind = reshape(l_self_ind,n_scl_pairs,size(l_self_ind,3)).';
l_ij_ji_ind = reshape(l_ij_ji_ind,n_gR,size(l_ij_ji_ind,3)).';
