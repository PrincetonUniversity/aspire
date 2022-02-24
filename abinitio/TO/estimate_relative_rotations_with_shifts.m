function est_rel_rots = estimate_relative_rotations_with_shifts(pf_norm, max_shift, shift_step, cache_file_name)

[n_r,~,n_images]= size(pf_norm);

shift_phases = calc_shift_phases(n_r, max_shift, shift_step);
n_shifts = size(shift_phases, 2);

g_shift_phases = gpuArray(single(shift_phases));

data = load(cache_file_name);
R = data.R;
n_candidates = size(R,3);
S_self = zeros(n_images, n_candidates);
est_rel_rots  = zeros(n_images, n_images);

n_pairs = nchoosek(n_images, 2);
ii_inds = zeros(n_pairs, 1);
jj_inds = zeros(n_pairs, 1);
clmats = cell(n_pairs, 1);

l_ij_ji_ind = data.l_ij_ji_ind;
l_self_ind = data.l_self_ind;
g_l_ij_ji_ind = gpuArray(single(l_ij_ji_ind));
g_l_self_ind = gpuArray(single(l_self_ind));

ind=1;
for ii = 1:n_images
    for jj = ii+1:n_images
        ii_inds(ind) = ii;
        jj_inds(ind) = jj;
        ind = ind+1;
    end
end

% self common lines
for p_i=1:n_images
    pf_norm_i = pf_norm(:,:,p_i);
    S_self(p_i,:) = gather(scls_correlation(n_shifts, g_shift_phases, pf_norm_i, g_l_self_ind));
end

% common lines
for ind = 1:n_pairs
    p_i = ii_inds(ind);
    p_j = jj_inds(ind);
    clmats{ind} = gather(max_correlation_pair_ind(n_shifts, g_shift_phases, pf_norm, p_i, p_j, S_self, n_candidates, g_l_ij_ji_ind));    
end

% store the rotations indices
for ind = 1:n_pairs
    p_i = ii_inds(ind);
    p_j = jj_inds(ind);
    [est_rel_rots(p_i,p_j),est_rel_rots(p_j,p_i)] = ind2sub([n_candidates n_candidates],clmats{ind});
end


function ind = max_correlation_pair_ind(n_shifts, g_shift_phases, pf_norm, p_i, p_j, S_self, n_candidates, l_ij_ji_ind)

Corrs_cls = zeros(size(l_ij_ji_ind),'gpuArray');

pf_norm_i = pf_norm(:,:,p_i);
g_pf_norm_i = gpuArray(single(pf_norm_i));
pf_norm_j = pf_norm(:,:,p_j);
g_pf_norm_j = gpuArray(single(pf_norm_j));

Sij = S_self(p_i, :)'*S_self(p_j, :);                   % self common lines correlation (n_candidates x n_candidates)
Sij(1:(n_candidates + 1):end) = 0;                      % Reset the diagonal

for ds=1:n_shifts  
    g_pf_norm_shifted = bsxfun(@times,g_pf_norm_j,g_shift_phases(:,ds));
    Corrs   = real(g_pf_norm_i'*g_pf_norm_shifted);
    Corrs_cls_tmp = Corrs(l_ij_ji_ind);
    Corrs_cls = max(Corrs_cls, Corrs_cls_tmp);
end

cl  = prod(Corrs_cls,2);
c  = reshape(cl, [n_candidates, n_candidates]);         % common lines correlation (n_candidates x n_candidates)

Sij = Sij.*c;                                           % element-wise multiplication

[~,ind] = max(Sij(:));                            % find the pair with the max correlation 


function res = scls_correlation(n_shifts, g_shift_phases, pf_norm_i, l_self_ind)

Corrs_scls = zeros(size(l_self_ind),'gpuArray');
g_pf_norm_i = gpuArray(single(pf_norm_i));

for ds=1:n_shifts
    g_pf_norm_shifted = bsxfun(@times,g_pf_norm_i,g_shift_phases(:,ds));
    Corrs_pi = real(g_pf_norm_i'*g_pf_norm_shifted);
    Corrs_scls_tmp = Corrs_pi(l_self_ind);
    Corrs_scls = max(Corrs_scls, Corrs_scls_tmp);
end
res = prod(Corrs_scls,2);
