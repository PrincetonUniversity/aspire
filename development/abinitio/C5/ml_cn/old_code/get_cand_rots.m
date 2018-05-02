function [cijs_inds,Ris_tilde,R_theta_ijs,n_theta] = get_cand_rots(n_symm,nPoints_sphere,n_theta,inplane_rot_res,is_use_gt_in_cands,refq)
% TODO: handle input variables

% TODO: have file name be an input variable and not hardcoded
file_cache_name  = sprintf('ml_cn_cached_inds_%d.mat',nPoints_sphere);
% file_cache_name  = sprintf('ml_c12_cached_inds_no_eq_%d.mat',nPoints_sphere);
% file_cache_name  = sprintf('ml_c12_cached_inds_%d_res%d.mat',nPoints_sphere,inplane_rot_res_cache);

if exist(file_cache_name, 'file') == 2
    log_message('loading line indeces cache');
    load(file_cache_name);
else
    % TODO: ask the user permission to calc and her permission to save to
    % file...
    log_message('generating line indeces cache');
    [Ris_tilde,R_theta_ijs] = generate_cand_rots(nPoints_sphere,inplane_rot_res,is_use_gt_in_cands,refq);
    cijs             = compute_cls_inds(Ris_tilde,R_theta_ijs,n_theta,is_save_inds_to_cache,file_cache_name);
end
    
n_theta_ijs = size(R_theta_ijs,3);

n_theta_ijs_to_keep = floor(n_theta_ijs/n_symm)*n_symm;
if n_theta_ijs_to_keep < n_theta_ijs
    cijs(:,:,n_theta_ijs_to_keep+1:end,:) = [];
    R_theta_ijs(:,:,n_theta_ijs_to_keep+1:end) = [];
end

cijs_inds = uint16(sub2ind([n_theta/2,n_theta],cijs(:,:,:,1),cijs(:,:,:,2)));