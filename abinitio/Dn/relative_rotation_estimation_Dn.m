function all_relative_rotations = relative_rotation_estimation_Dn(...
    symmetry_degree, pf_images, gR, cands)
%RELATIVE_ROTATION_ESTIMATION_DN Generation of Dihedral group elements
%   This function generates the elements of the dihedral symmetry-group Dn.
%   Input:
%       symmetry_degree - an integer >= 3.
%       pf_images - Polar Fourier-Transformed projection images. The i-th
%                   image is given by pf_images(:,:,i).
%       gR - Group elements of Dn, gR(:,:,i) is the i-th group-element.
%       cands - For debugging only.
%
%   Output:
%       all_relative_rotations - All relative-rotations (R^T_i*g*R_j)
%                                for all g in Dn and i<j indices. 
%                                The k-th relative-rotation of the p-th
%                                pair (in linear indices) is given by
%                                all_relative_rotations(p,:,:,k).
%       
%   Written by Elad Eatah April 2021. 

[n_r, n_theta,~] = size(pf_images);

%% Candidates Cache Fetching
% Looking for an existing cache of candidates matrices and all 
% their common-lines (including self common-lines). If no such cache is
% available, for our specific symmetry degree (n), then such cache is
% generated. For debugging, this cache is generated from the matrices
% 'cands'.
folder_name = 'Candidates sets';
cache_file_name = fullfile(folder_name, 'D' + string(symmetry_degree) + ...
    '_symmetry*.mat');
cache_struct = dir(cache_file_name);
if isempty(cache_struct)
    log_message('Candidates cache for D' + string(symmetry_degree) + ...
        ' not found')
    log_message('Generating candidates cache');
    cache_file_name = cryo_Dn_create_cache(symmetry_degree, ...
        gR(:,:,1+symmetry_degree), cands, folder_name);
else
    cache_file_name = fullfile(cache_struct.folder, cache_struct.name);
end
cache_vars = load(cache_file_name);
candidates = cache_vars.R;
[candidates_cijs, candidates_cjis]  = ind2sub([n_theta, n_theta],...
    cache_vars.l_ij_ji_ind);
[candidates_ciis, candidates_cigs] = ind2sub([n_theta, n_theta],...
    cache_vars.l_self_ind);

%% Configuration
max_shift = 1;
shift_step = 0.5;
%precompile the shift phases and normalizing all input images.
shift_phases = calc_shift_phases(n_r,max_shift,shift_step);
normalized_pf_images = normalize_images(pf_images);
rotations_num = size(pf_images, 3);  % N - number of rotations to estimate.
candidates_num = size(candidates, 3);  % L - number of candidates in grid.

rotation_pairs_count = nchoosek(rotations_num, 2); % Pairs count {i<j}
candidate_pairs_count = candidates_num ^ 2; % All candidates pairs.
self_common_lines_indices = 1:(candidates_num + 1):candidate_pairs_count;
common_lines_scores = ones(rotation_pairs_count, candidate_pairs_count);
self_common_lines_scores = NaN(rotations_num, candidates_num);

group_order = size(gR, 3);  % 2n for the group Dn.

%% Utilizing self common-lines.
% Estimating the correlation of self common-lines (ro_{ii} in the
% algorithm).
for rotation_ind = 1:rotations_num
    npf_image_i = normalized_pf_images(:,:,rotation_ind);
         
    % corrs is candidates_num X (common_lines_count - 1), since for self
    % common-lines one does not use the identity element in the group.
    corrs = common_lines_correlation(...
         candidates_ciis, candidates_cigs, npf_image_i, npf_image_i, ...
         shift_phases);
    self_common_lines_scores(rotation_ind, :) = prod(real(corrs), 2);
end

%% Utilizing all other common-lines.
% Estimating using all common-lines and the previously-computed
% correlation of self common-lines (pi_{ij} in the algorithm).
for rotation_ind = 1:rotations_num    
    npf_image_i = normalized_pf_images(:,:,rotation_ind);
    
    for other_rotation_ind = rotation_ind+1:rotations_num
        pair_ind = uppertri_ijtoind_vec(rotation_ind, ...
            other_rotation_ind, rotations_num);
        npf_image_j = normalized_pf_images(:,:,other_rotation_ind);
        
        corrs = common_lines_correlation(...
            candidates_cijs, candidates_cjis, npf_image_i, npf_image_j, ...
            shift_phases);  % corrs is pairs_count X common-lines_count.
        common_lines_scores(pair_ind, :) = prod(real(corrs), 2);
        
        for candidate_pair = 1:candidate_pairs_count
            [candidate_l, candidate_k] = ind2sub(...
                [candidates_num, candidates_num], candidate_pair);
            common_lines_scores(pair_ind, candidate_pair) = ...
                common_lines_scores(pair_ind, candidate_pair) * ...
                self_common_lines_scores(rotation_ind, candidate_l) * ...
                self_common_lines_scores(other_rotation_ind, candidate_k);
        end        
    end
end

% Assuming all rotations are unique, we enforce that
% self relative rotations can NOT be chosen for relative rotations
% for i<j.
common_lines_scores(:, self_common_lines_indices) = -1;

%% Maximum-Likelihood Estimation using all common-lines.
% Picking the candidates which attain the maximum total correlation score
[~, max_indices] = max(common_lines_scores, [], 2);
[max_indices_l, max_indices_k] = ind2sub(...
    [candidates_num, candidates_num], max_indices);
estimated_Ri = candidates(:, :, max_indices_l);
estimated_Rj = candidates(:, :, max_indices_k);

% Constructing all relative-rotations for each candidates-pair.
all_relative_rotations = NaN(rotation_pairs_count, 3, 3, group_order);
for pair_ind = 1:rotation_pairs_count
    all_relative_rotations(pair_ind, :, :, :) = multiprod(multiprod(...
        estimated_Ri(:, :, pair_ind).', gR), estimated_Rj(:, :, pair_ind));
end
end


function corrs = common_lines_correlation(cijs, cjis, ...
    image_i, image_j, shift_phases)
%COMMON_LINES_CORRELATION Computes the common-lines correlation score
%between two input images and estimated common-lines.
%   Input:       
%       cijs - common-lines in the coordinates system of the 
%              first image's projection-plane.
%       cjis - common-lines in the coordinates system of the 
%              second image's projection-plane.
%       image_i - First image, after Polar-Fourier Transform and
%                 normalization.
%       image_j - Second image, after Polar-Fourier Transform and
%                 normalization.
%       shift-phases - A list of shifts for correlation computation.
%
%   Output:
%       corrs - 2D array of correlations-scores.
%       
%   Written by Elad Eatah April 2021. 

    % Creating all shifts of image_{i}.
    [n_r, n_theta] = size(image_i);
    nshifts = size(shift_phases, 2);
    common_lines_count = size(cijs, 2);
    npf_i_shifted = zeros(n_r,n_theta, nshifts);
    for s=1:nshifts
        npf_i_shifted(:,:,s) = bsxfun(@times,image_i,shift_phases(:,s));
    end
    
    % Normalizing all shifted copies of image_{i}.
    npf_i_shifted = normalize_images(npf_i_shifted);
    Npf_i_shifted = gpuArray(single(npf_i_shifted)); %single(npf_i_shifted);
    Npf_j = gpuArray(single(image_j)); %single(image_j);
    
    candidates_num = size(cijs, 1);
    
    % Computing correlation scores.
    % cross-correaltion - so we want to conjugate
    co = bsxfun(@times, Npf_i_shifted(:,cijs,:), conj(Npf_j(:,cjis))); 
    Corrs = sum(co);
    corrs = gather(Corrs);        
    corrs = reshape(corrs, candidates_num, common_lines_count, nshifts);
    
    if nshifts > 1
        % the 1d shifts depend not only on the 2d shifts of each and every
        % image, but also on their common-lines. That is, for a given image, the 1d
        % shift might be different for different lines !
        % In particular, this means that each one of the 4 lines we are considering in a given image
        % may correspond to a DIFFERENT 1d shift. Therefore each line has the dof, so to speak, of choosing its own shift.
        % All in all, this means that we max out the shift BEFORE taking
        % thet mean score over all n lines.
        corrs = max(real(corrs),[],3);
    end
end

function normalized_images = normalize_images(pf_images)
%NORMALIZE_IMAGES Normalizes all polar-images per-ray, such that
% each ray (i.e row) has L2-norm 1.
%   Input:       
%       pf_images - Polar Fourier-Transformed projection images. The i-th
%                   image is given by pf_images(:,:,i).
%
%   Output:
%       normalized_images - 3D array of normalized images, the i-th image
%                           is given by normalized_images(:,:,i).
%       
%   Written by Elad Eatah April 2021. 

    % ignoring dc term
    npf_tmp = pf_images;
    npf_tmp(1,:,:) = 0;
    normalized_images = bsxfun(@rdivide,npf_tmp,...
        sqrt(sum((abs(npf_tmp)).^2)));
end
