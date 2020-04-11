function estimatedVol = reconstruct_vol_D2(noisy_projs,npf,rot_alligned,params)

[noisy_projs_replicated,npf_replicated]  = replicate_projs(noisy_projs,npf);
if ~exist('results','dir')
    mkdir('results');
end
rot_alligned_g_duplicated = g_duplicate_rots(rot_alligned);

estimatedVol = reconstruct3d(npf_replicated,noisy_projs_replicated,rot_alligned_g_duplicated,params);
WriteMRC(estimatedVol,1,'results/1.mrc');

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

global g;
[h,w,nRots] = size(rots);
assert(h == 3 && w == 3);

rots_g_duplicated = zeros(3,3,4*nRots);
gx=diag([1,-1,-1]);
gy=diag([-1,1,-1]);
gz=diag([-1,-1,1]);
for k=1:nRots
    rot_k = rots(:,:,k);
    rots_g_duplicated(:,:,4*(k-1)+1) = rot_k;
    rots_g_duplicated(:,:,4*(k-1)+2) = gx*rot_k;
    rots_g_duplicated(:,:,4*(k-1)+3) = gy*rot_k;
    rots_g_duplicated(:,:,4*(k-1)+4) = gz*rot_k;
    
%     for s=1:4
%         rots_g_duplicated(:,:,4*(k-1)+s) = g^(s-1)*rot_k;
%     end
end

end


function vol = reconstruct3d(npf,projs,rots,params)

% %TODO: Gabi remove me
% rots_gt = zeros(3,3,params.K);
% for i=1:params.K
%     rots_gt(:,:,i) = q_to_rot(params.refq(:,i))';
% end
% [est_shifts,~] = cryo_estimate_shifts(npf,rots,ceil(2*sqrt(2)*params.max_shift),...
%                                       params.shift_step,10000,[],1);
% [est_shifts,~] = cryo_estimate_shifts(npf,rots,nshifts,params.shift_step,10000,[],0);
[est_shifts,~] = cryo_estimate_shifts(npf,rots,params.max_shift,...
                                      params.shift_step,10000,[],1);

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
