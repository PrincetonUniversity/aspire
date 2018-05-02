function estimatedVol = reconstruct_vol(noisy_projs,npf,rot_alligned,max_shift,shift_step)

if ~exist('shift_step','var')
    shift_step = 0.5;
end

if ~exist('max_shift','var')
    max_shift = 15;
end

log_message('reconstructing volume');

[noisy_projs_replicated,npf_replicated]  = replicate_projs(noisy_projs,npf);
rot_alligned_g_duplicated = g_duplicate_rots(rot_alligned);

estimatedVol = reconstruct3d(npf_replicated,noisy_projs_replicated,rot_alligned_g_duplicated,max_shift,shift_step);

end


function [projs_replicated,npf_replicated] = replicate_projs(projs,npf)

[h1,w1,nProjs1] = size(projs);
[h2,w2,nProjs2] = size(npf);

assert(nProjs1 == nProjs2);

projs_replicated = repmat(projs,[1 4 1]);
projs_replicated = reshape(projs_replicated,[h1,w1,4*nProjs1]);

npf_replicated = repmat(npf,[1 4 1]);
npf_replicated = reshape(npf_replicated,[h2,w2,4*nProjs2]);


end


function rots_g_duplicated = g_duplicate_rots(rots)

g = [0 -1 0; ...
     1  0 0; ...
     0  0 1]; % rotation matrix of 90 degress about z-axis

[h,w,nRots] = size(rots);
assert(h == 3 && w == 3);

rots_g_duplicated = zeros(3,3,4*nRots);
for k=1:nRots
    rot_k = rots(:,:,k);
    for s=1:4
        rots_g_duplicated(:,:,4*(k-1)+s) = g^(s-1)*rot_k;
    end
end

end


function vol = reconstruct3d(npf,projs,rots,max_shift,shift_step)

[est_shifts,~] = cryo_estimate_shifts(npf,rots,ceil(2*sqrt(2)*max_shift),...
                                      shift_step,10000,[],1);
% [est_shifts,~] = cryo_estimate_shifts(npf,rots,nshifts,params.shift_step,10000,[],0);

% if ~params.real_data && params.debug && isfield(params,'ref_shifts')
%     [est_shifts,~] = cryo_estimate_shifts(npf,rots_inv,nshifts,...
%         params.shift_step,10000,params.ref_shifts.',0);
% else
%     [est_shifts,~] = cryo_estimate_shifts(npf,rots_inv,nshifts,params.shift_step,10000,[],0);
% end

n = size(projs,1);
[ vol, ~, ~ ,~, ~, ~] = recon3d_firm( projs,...
    rots,-est_shifts, 1e-6, 100, zeros(n,n,n));

% [ vol, ~, ~ ,~, ~, ~] = recon3d_firm( projs,...
%     rots,[], 1e-6, 100, zeros(n,n,n));


vol = real(vol);

% [vol,~,~,~,~,~] = recon3d_firm(projs,...
%     rots,[], 1e-6, 30, zeros(89,89,89));

% [vol,~,~,~,~,~] = recon3d_firm(projs,rots,[]);
% 
% assert(norm(imag(vol(:)))/norm(vol(:))<1.0e-3);
% vol = real(vol
end