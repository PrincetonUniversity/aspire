function [images,refq,ref_shifts,clean_images,vol] = generate_cn_images(n_symm,nImages,SNR,projSize,cn_type,max_shift,shift_step)
%
% Generates a set of projection images of a C4 volume
% 
% Input parameters:
%   nImages       Number of images to generate
%   SNR           Signal to noise ratio
%   projSize      Size of each projection (e.g., 65 would generate a 65x65 images)
%   c5_type       (Optional) Type of c4 molecule. May be
%                 'GAUSSIAN','SYNTHETIC','IP3','TRPV1'. Default is GAUSSIAN.
%   max_shift     The maximum spatial shift each image undergoes
%   shift_step    (Optional) Resolution used to generate shifts.
%                  shift_step=1 allows for all integer shifts from
%                  -max_shift to max_shift (in both x and y directions).
%                  shift_step=2 generates shifts between -max_shift and
%                  max_shift with steps of 2 pixels, that is, shifts of 0
%                  pixels, 2 pixels, 4 pixels, and so on. Default
%                  shift_step is 0. 
%
% Output parameters:
%   images           A 3-dim array of images of size projSize x projSize x nImages
%   refq             A array of size 4xnImages. The i-th column is a
%                    quaternion representing the rotation matrix of the
%                    i-th image.
%   ref_shifts       A two column table with the 2D shift introduced to
%                    each projections.
              

if ~exist('shift_step','var')
    shift_step = 0.1;
end

if ~exist('max_shift','var')
    max_shift = 15;
end

log_message('generating images of class C%d',n_symm);

refq   = qrand(nImages);
rots = q_to_rot(refq);

if strcmp(cn_type,'80S')
%     vol_init_mat = cryo_fetch_emdID(2660);
    vol_init_mat = '/tmp/tp8641e8a4_c426_4a2d_8028_6e2d046b7cb0/pub/databases/emdb/structures/EMD-2660/map/emd_2660.map';
    vol_init = ReadMRC(vol_init_mat);
    vol_init = cryo_downsample(vol_init,89);
    vol_lp = GaussFilt(vol_init,0.15);
    masked_vol = cryo_mask_volume(vol_lp,25,5);
    vol = make_vol_cn_old(n_symm,masked_vol);
    vol = cryo_mask_volume(vol,25,5);
    assertVolumeIsCn(vol,n_symm);
elseif strcmp(cn_type,'80S_DUPLICATED')
%     vol_mat = cryo_fetch_emdID(2660);
    vol_mat = '/tmp/tp8641e8a4_c426_4a2d_8028_6e2d046b7cb0/pub/databases/emdb/structures/EMD-2660/map/emd_2660.map';
    vol = ReadMRC(vol_mat);
    vol = make_vol_cn(vol,n_symm);
    assertVolumeIsCn(vol,n_symm);
elseif strcmp(cn_type,'C1')
    vol = cryo_gaussian_phantom_3d('C1_params',projSize,1);
    vol = make_vol_cn_old(vol);
    assertVolumeIsCn(vol,n_symm);
    WriteMRC(vol,1,sprintf('ml_c%d_init.mrc',cn_type));
elseif strcmp(cn_type,'C1_DUPLICATED')
    vol = cryo_gaussian_phantom_3d('C1_params',projSize,1);
    vol = make_vol_cn(vol,n_symm);
    assertVolumeIsCn(vol,n_symm);
else
    error('no such type %s',cn_type);
end

log_message('#images = %d', nImages);

% projs = cryo_project_gaussian('C5_params',projSize,1,refq);
projs = cryo_project(vol,rots);
projs = permute(projs,[2,1,3]);

log_message('max shifts=%d, shift step=%7.5f',max_shift,shift_step);
[projs,ref_shifts] = cryo_addshifts(projs,[], max_shift,shift_step);
% [projs,params.ref_shifts] = cryo_addshifts(projs,2*ones(params.nImages,2));

log_message('SNR=%7.5f',SNR);
images = cryo_addnoise(projs,SNR,'gaussian');

clean_images = projs;

%     [~,noisy_projs,~,params.refq] = cryo_gen_projections(params.projSize,params.nImages,params.SNR,params.max_shift,1,...
%         [],params.refq,'single');
%     [~,noisy_projs,~,params.refq] = cryo_gen_projections_vol(symmVol,params.projSize,params.nImages,params.SNR,...
%         params.max_shift,1,...
%         [],params.refq,'single');
%     spnoise_pft = cryo_noise_estimation_pfs(noisy_projs,params.n_r,params.n_theta);
%     spnoise_pft(spnoise_pft<max(spnoise_pft/10))=1;
%     noise_cov = diag(spnoise_pft);

end

function vol_out = make_vol_cn(vol_in,n_symm)

inn = size(vol_in,1);

vol_tmp = zeros(24+inn,24+inn,24+inn);
% 
% bttm = -floor(inn/2);
% top  = floor(inn/2);
vol_tmp(1:inn,1:inn,1:inn) = vol_in;

vol_out = vol_tmp;
for i=1:n_symm-1
    vol_out = vol_out + fastrotate3z(vol_tmp,i*360/n_symm);
end
vol_out = vol_out/n_symm;

end



function vol_out = make_vol_cn_old(n_symm,vol_in)

vol_out = vol_in;
for i=1:n_symm-1
    vol_out = vol_out + fastrotate3z(vol_in,i*360/n_symm);
end
vol_out = vol_out/n_symm;

end


function err = assertVolumeIsCn(vol,n_symm)
   
errs = zeros(1,n_symm);
for i=1:n_symm-1
    vol_rot = fastrotate3z(vol,360*i/n_symm);
    errs(i) = norm(vol(:)-vol_rot(:))/norm(vol(:));
end
err = max(errs);
log_message('deviation of volume from c%d symmetry is %.4f',n_symm,err);
end


% function refq = removeEquators(refq)
% 
% log_message('removing equator images');
% nImages = size(refq,2);
% is_eq = zeros(1,nImages);
% for i=1:nImages
%     rot = q_to_rot(refq(:,i)).';
%     if abs(acosd(rot(3,3))-90) < 8
%         is_eq(i) = 1;
%     end
% end
% 
% refq(:,find(is_eq)) = [];
% 
% end