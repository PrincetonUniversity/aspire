function [images,refq,ref_shifts] = generate_c2_images(nImages,SNR,projSize,c2_type,max_shift,shift_step)
%
% Generates a set of projection images of a C4 volume
% 
% Input parameters:
%   nImages       Number of images to generate
%   SNR           Signal to noise ratio
%   projSize      Size of each projection (e.g., 65 would generate a 65x65 images)
%   c4_type       (Optional) Type of c4 molecule. May be
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

if ~exist('c2_type','var')
    c2_type = 'GAUSSIAN';
end

refq   = qrand(nImages);

log_message('loading simulated %s volume',c2_type)
if strcmp(c2_type,'GAUSSIAN') 
    vol = cryo_gaussian_phantom_3d('C2_params',projSize,1);   
elseif strcmp(c2_type,'SYNTHETIC')
    
    load cleanrib;
    volref = real(volref);
    vol = volref + fastrotate(volref,180);
    clear volref;
else
    error(['no such c4_type : ' c2_type]);
end

% if params_simul.debug
%     %     assertVolumeCentered(symmVol);
%     %     assertVolumeIsC4(symmVol);
% end

log_message('#images = %d', nImages);

rots = q_to_rot(refq);
projs = cryo_project(vol,rots);
projs = permute(projs,[2,1,3]);

log_message('max shifts=%d, shift step=%7.5f',max_shift,shift_step);
[projs,ref_shifts] = cryo_addshifts(projs,[], max_shift,shift_step);
% [projs,params.ref_shifts] = cryo_addshifts(projs,2*ones(params.nImages,2));

log_message('SNR=%7.5f',SNR);
images = cryo_addnoise(projs,SNR,'gaussian');

%     [~,noisy_projs,~,params.refq] = cryo_gen_projections(params.projSize,params.nImages,params.SNR,params.max_shift,1,...
%         [],params.refq,'single');
%     [~,noisy_projs,~,params.refq] = cryo_gen_projections_vol(symmVol,params.projSize,params.nImages,params.SNR,...
%         params.max_shift,1,...
%         [],params.refq,'single');
%     spnoise_pft = cryo_noise_estimation_pfs(noisy_projs,params.n_r,params.n_theta);
%     spnoise_pft(spnoise_pft<max(spnoise_pft/10))=1;
%     noise_cov = diag(spnoise_pft);

end