function cryo_abinitio_Cn_ml(n_symm,instack,outvol,outparams,...
    nImages,n_r,max_shift,shift_step)

nPoints_sphere  = 1000;
inplane_rot_res = 1;

initstate; 
open_log(0)

is_handle_equators = false;
if is_handle_equators == true
    assert(mod(n_symm,2) == 0);
end

if ~exist('n_r','var')
    % defer exact value once the size of images is known
    n_r_given = false;
else
    n_r_given = true;
end

if ~exist('max_shift','var')
    % defer maximum shift once the size of images is known
    max_shift_given = false;
else
    max_shift_given = true;
end

if ~exist('shift_step','var')
    shift_step = 0.5;
end

% Load projections
projs = ReadMRC(instack);
% projs = projs(:,:,floor(linspace(1,size(projs,3),nImages)));
% log_message('REMOVE THIS LINE!!!');
% projs = projs(:,:,1:5000);
% projs = projs(:,:,floor(linspace(1,size(projs,3),nImages)));

% im_inds = 1:nImages;

% projs = projs(:,:,1:5000);
im_inds = floor(linspace(1,size(projs,3),nImages));

% run this line only once in order to generate the candidate rots
% [cijs_inds,Ris_tilde,R_theta_ijs] = get_cand_rots(n_symm,nPoints_sphere,n_theta,inplane_rot_res,is_use_gt_in_cands,refq);
[cijs_inds,Ris_tilde,R_theta_ijs,n_theta] = get_cand_rots(n_symm,nPoints_sphere);
ciis = compute_scls_inds(Ris_tilde,n_symm,n_theta);


projs = projs(:,:,im_inds);
save(outparams,'im_inds');

nImages = size(projs,3);

log_message('symmetry class is C%d',n_symm);
log_message('projections loaded. Using %d projections of size %d x %d',nImages,size(projs,1),size(projs,2));
if size(projs,1)~=size(projs,2)
    error('Input images must be square');
end

figure; viewstack(projs,5,5);

mask_radius = round(size(projs,1)*0.4);
% mask_radius = 60;
log_message('Masking projections. Masking radius is %d pixels',mask_radius);
masked_projs = mask_fuzzy(projs,mask_radius);
save(outparams,'mask_radius','-append');
% 
% log_message('REMOVE THIS!!!!!!!');
% masked_projs = projs;
    
if ~n_r_given
    n_r = ceil(size(masked_projs,1)*0.5);
%     log_message('REMOVE THIS!!!!!!!');
%     n_r = 45;
    save(outparams,'n_r','-append');
end
[npf,~] = cryo_pft(masked_projs,n_r,n_theta,'double');

npf = gaussian_filter_imgs(npf);

if ~max_shift_given
    max_shift = ceil(size(projs,1)*0.15); % max_shift is 15% of the image size
%     max_shift = 3;
    save(outparams,'max_shift','-append');
end

[vijs,viis,max_corrs_stats] = compute_third_row_outer_prod_both2(npf,ciis,cijs_inds,Ris_tilde,R_theta_ijs,...
    n_symm,max_shift,shift_step,is_handle_equators);
% [vijs,viis,max_corrs_stats] = compute_third_row_outer_prod_both(npf,ciis,cijs,Ris_tilde,R_theta_ijs,max_shift,shift_step,is_handle_equators);

vijs = mat2flat(vijs,nImages);
[vijs,viis,sign_ij_J,sign_ii_J] = global_sync_J(vijs,viis);
% 
vis  = estimate_third_rows_ml(vijs,viis);
rots = estimate_inplane_rotations(npf,vis,n_symm,inplane_rot_res,max_shift,shift_step);

% [rot_alligned,err_in_degrees,mse] = analyze_results_ml(rots,n_theta,refq);


save(outparams,'rots','-append');

estimatedVol = reconstruct_ml_cn(masked_projs,rots,n_symm,n_r,n_theta,max_shift,shift_step);
WriteMRC(estimatedVol,1,outvol);


post_process(outvol,n_symm);

end