function cache_filename = cryo_Dn_create_cache(symmetry_degree,...
    order_2_gen, cands, files_path)

% Creates a mat file containing:
%       R - candidate rotations set. 
%       n_theta
%       l_ij_ji_ind - set of common lines in linear indices.
%       l_self_ind - set of self common lines in linear indices.
%
%   Input:
%       symmetry_degree - an integer >= 3.
%
%   Output:
%       cache_filename - a mat file name containing all required data.
%       
%   Written by Elad Eatah March 2021. 

n_theta = 360;
inplane_rot_degree = deg2rad(0.5);
n_scl_pairs = 2 * symmetry_degree - 1;  % No. of self common-lines.
n_gR = 2 * symmetry_degree;  % Group order
viewing_angle = -0.9;  % Cosine of the viewing angle
resolution = 30;
rot_axis = rotm2axang(order_2_gen);
rot_axis = rot_axis(1:3);

%[~, ~, ~, resolution, viewing_angle, inplane_rot_degree] = cryo_TO_configuration();
gR = cryo_Dn_group_elements(symmetry_degree, order_2_gen);

log_message('Generating candidate rotations set');
R = cands; %TO_candidate_rotations_set(gR, resolution, viewing_angle, inplane_rot_degree);
n_candidates = size(R,3);

log_message('Computing common lines and self common lines indices sets');
[l_self_ind, l_ij_ji_ind] = Dn_cl_scl_inds(gR, n_gR, n_scl_pairs, n_theta, R);

cache_filename = sprintf('D%d_symmetry_%d_candidates_cache_test.mat',symmetry_degree, n_candidates);
cache_filename = fullfile(files_path, cache_filename);
log_message('Cache file name: %s', cache_filename);
save(cache_filename,'l_self_ind','l_ij_ji_ind','-v7.3','-nocompression');
save(cache_filename,'R','n_theta','rot_axis','-append');    
end

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
end

function [l_self_ind, l_ij_ji_ind] = Dn_cl_scl_inds(gR, n_gR, n_scl_pairs, n_theta, R)

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
%                       n_candidates_pairs x n_gR
%       l_ij_ji_ind -   set of common lines in linear indices:
%                       n_candidates x n_scl_pairs
%       
%   Written by Adi Shasha January 2021. 
%   Modified by Elad Eatah for Dn April 2021.

n_R = size(R, 3);

% array of linear indices of common lines and self common lines
l_self  = ones(2,n_scl_pairs,n_R);
l_ij_ji = ones(2,n_gR,n_R*n_R);

[all_cls1, all_cls2] = find_all_common_lines(R, gR, n_theta);

    
for i=1:n_R 
    % Compute self common lines of Ri.
    l_self(1,:,i) = permute(all_cls1(i, i, 2:n_gR), [1, 3, 2]);    
    l_self(2,:,i) = permute(all_cls2(i, i, 2:n_gR), [1, 3, 2]);
    a = NaN('like', l_self);
%     for k=2:n_scl_pairs
%         %[l_self(1,k,i),l_self(2,k,i)] = commonline_R(R(:,:,i),R(:,:,i)*gR(:,:,k),n_theta);
%         [a(1,k,i),~] = commonline_R(R(:,:,i),R(:,:,i)*gR(:,:,k),n_theta);
%     end
    for j=1:n_R 
        if i == j
            continue
        end        
        % Compute the common lines induced by rotation matrices Ri and Rj.
        ind = sub2ind([n_R,n_R],i,j);
        l_ij_ji(1,:,ind) = permute(all_cls1(i, j, :), [1, 3, 2]); 
        l_ij_ji(2,:,ind) = permute(all_cls2(i, j, :), [1, 3, 2]);
%         ind = sub2ind([n_R,n_R],i,j);
%          for k=1:n_gR
%              [l_ij_ji(1,k,ind),l_ij_ji(2,k,ind)] = commonline_R(R(:,:,i),R(:,:,j)*gR(:,:,k),n_theta);
%          end
    end
end

%l_self  = l_self + 1;
%l_ij_ji = l_ij_ji + 1;

% convert to linear indices
l_self_ind  = sub2ind([n_theta,n_theta],l_self(1,:,:),l_self(2,:,:));
l_ij_ji_ind = sub2ind([n_theta,n_theta],l_ij_ji(1,:,:),l_ij_ji(2,:,:));

% reshape
l_self_ind = reshape(l_self_ind,n_scl_pairs,size(l_self_ind,3)).';
l_ij_ji_ind = reshape(l_ij_ji_ind,n_gR,size(l_ij_ji_ind,3)).';
end


function [cands_scl_cijs, cands_scl_cjis] = find_all_common_lines(...
    cands, gR, n_theta)
    candidates_num = size(cands,3);
    group_order = size(gR, 3);
    
    cands_scl_cijs = NaN(candidates_num, candidates_num, group_order);
    cands_scl_cjis = NaN(candidates_num, candidates_num, group_order);
    
    for cand_ind = 1:candidates_num
        % Separate computation for self common_lines (only 2n-1 lines)
        [c1, c2] = get_common_lines(cands(:,:,cand_ind), ...
                cands(:,:,cand_ind), gR(:,:,2:end), n_theta);
            cands_scl_cijs(cand_ind, cand_ind, 2:end) = c1;
            cands_scl_cjis(cand_ind, cand_ind, 2:end) = c2;
        
        for other_cand_ind = [1:cand_ind-1, cand_ind+1:candidates_num]
            [c1, c2] = get_common_lines(cands(:,:,cand_ind), ...
                cands(:,:,other_cand_ind), gR, n_theta);
            cands_scl_cijs(cand_ind, other_cand_ind, :) = c1;
            cands_scl_cjis(cand_ind, other_cand_ind, :) = c2;
        end
    end
end

function [cijs, cjis] = get_common_lines(Ri, Rj, group_elements, n_theta)
    group_elements_count = size(group_elements,3);
    Us  = multiprod(multiprod(Ri',group_elements), Rj);
    Us = reshape(Us, 9, group_elements_count);
    % extract a common-line index for each possible theta_ij
    c_1 = [-Us(8,:) ;  Us(7,:)]; %[-U(2,3)  U(1,3)];
    c_2 = [ Us(6,:) ; -Us(3,:)]; %[ U(3,2) -U(3,1)];

    cijs = clAngles2Ind(c_1, n_theta);
    cjis = clAngles2Ind(c_2, n_theta);
end

