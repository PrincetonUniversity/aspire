function vijs_out = HandednessSynchronizationDn(relative_rotations, ...
    rotations_num, symmetry_degree)
%HandednessSynchronizationDn Estimates the outer-products of the
%   desired-rotations' third-rows. Either all output products of spurios-J in
%   them or none have such J.
%
%   Input:
%       relative_rotations - Estimated relative-rotations (calculated in
%                            the former step of the algorithm).
%       rotations_num - The number of orientations to estimate (N).
%       symmetry_degree - an integer >= 3 (n).
%
%   Output:
%       vijs_out - All such outer products as 3X3X(N choose 2) matrix.
%       
%   Written by Elad Eatah March 2021.

% Pick all estimations of the requested outer-products,as an (N choose
% 2)X3X3X2 array (each product is estimated twice).
vijs_cands = fetch_third_rows_outer_products(relative_rotations, ...
    rotations_num, symmetry_degree);

% Select the closest approximations to rank-1 matrices
vijs = pick_best_candidates(vijs_cands);  % 3X3X(n choose 2).

% Form the handedness-graph and perform the synchronization.
[vijs_out,sign_ij_J] = global_sync_J_Dn(vijs, rotations_num);
end


function vijs_cands = fetch_third_rows_outer_products(...
    relative_rotations, rotations_num, symmetry_degree)
%FETCH_THIRD_ROWS_OUTER_PRODUCT Estimates the outer-products of the
%   desired-rotations' third-rows. Each product is estimated twice, by
%   averaging on Cn's relative-rotations and by averaging on Dn-Cn's
%   relative-rotations (with minus).
%
%   Input:
%       relative_rotations - Estimated relative-rotations (calculated in
%                            the former step of the algorithm).
%       rotations_num - The number of orientations to estimate (N).
%       symmetry_degree - an integer >= 3 (n).
%
%   Output:
%       vijs_cands - All such outer products estimations as an (N choose
%                    2)X3X3X2 array.
%       
%   Written by Elad Eatah March 2021.

pairs_count = nchoosek(rotations_num, 2);
vijs_cands = NaN(pairs_count,3,3,2);

for i=1:rotations_num
    ri_indices = uppertri_ijtoind_vec(i, i+1:rotations_num, rotations_num);
    for j=i+1:rotations_num
        pair_ind = ri_indices(j - i);
        rijs = squeeze(relative_rotations(pair_ind,:,:,:));
        vijs_cands(pair_ind,:,:,1) = mean(rijs(:,:,1:symmetry_degree), 3);
        vijs_cands(pair_ind,:,:,2) = -mean(rijs(:,:,symmetry_degree+1:2*symmetry_degree), 3);
    end
end

end


function vijs = pick_best_candidates(vijs_cands)
%PICK_BEST_CANDIDATES Selects the candidate, for each outer-product, which
%   is the closest to a rank-1 matrix.
%
%   Input:
%       vijs_cands - A (N choose 2)X3X3X2 array of outer-products'
%                    estimations.
%
%   Output:
%       vijs - Approximated rank-1 outer products as 3X3X(N choose 2) matrix.
%       
%   Written by Elad Eatah March 2021.

elements_count = size(vijs_cands,1);
candidates_count = size(vijs_cands,4);
vijs = NaN(3,3,elements_count);

for i = 1:elements_count
    cands = squeeze(vijs_cands(i,:,:,:));
    cands_scores = NaN(1, candidates_count);
    for j=1:candidates_count
        this_cand = squeeze(cands(:,:,j));
        sigmas = svds(this_cand);
        sigmas(1) = sigmas(1) - 1;
        cands_scores(j) = norm(sigmas);
    end
    [~,best_cand_index] = min(cands_scores); 
    vijs(:,:,i) = cands(:,:,best_cand_index);
end
end
