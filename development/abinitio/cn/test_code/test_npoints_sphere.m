function [detec_rate_err,results_mse] = test_npoints_sphere()

nPoints_sphere_arr = [2000];
nImages_arr = [200];

initstate; 
open_log(0)

results_mse    = zeros(numel(nImages_arr),numel(nPoints_sphere_arr));
detec_rate_err = zeros(numel(nImages_arr),numel(nPoints_sphere_arr));
for i=1:numel(nImages_arr)
    
   nImages = nImages_arr(i);
   
   for j=1:numel(nPoints_sphere_arr)
       
        nPoints_sphere = nPoints_sphere_arr(j);
       
        [mse,detec_rate] = do_test_npoints_sphere(nImages,nPoints_sphere);
        
        results_mse(i,j) = mse;
        detec_rate_err(i,j) = detec_rate;
        
        log_message('nImages=%d,nPoints_sphere=%d: mse=%.2f, detec_rate=%.2f%%', ...
                     nImages,   nPoints_sphere, mse, detec_rate*100);
   end
end

end

function [mse_vij,detec_rate] = do_test_npoints_sphere(nImages,nPoints_sphere)

is_use_gt_in_cands = false;
is_save_inds_to_cache = true;
n_r     = 65;  
n_theta = 360;
snr = 1000000000;
inplane_rot_res = 1;
c5_type = 'C1';
is_removeEquators = true;
max_shift  = 0;
shift_step = 0.5;

[projs,refq] = generate_c5_images(nImages,snr,65,c5_type,is_removeEquators,max_shift,shift_step);
nImages = size(refq,2);

[npf,~]      = cryo_pft(projs,n_r,n_theta,'double');

file_cache_name  = sprintf('ml_c5_cached_inds_%d.mat',nPoints_sphere);

if is_use_gt_in_cands || is_save_inds_to_cache
    [Ris_tilde,R_theta_ijs] = generate_cand_rots(nPoints_sphere,inplane_rot_res,is_use_gt_in_cands,refq);
    [ciis,cijs]             = compute_cls_inds(Ris_tilde,R_theta_ijs,n_theta,is_save_inds_to_cache,file_cache_name);
else
    load(file_cache_name);
end

[~,~,~,mse_vij,detec_rate] = compute_third_row_outer_prod(npf,ciis,cijs,Ris_tilde,R_theta_ijs,max_shift,shift_step,refq);


end