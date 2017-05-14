function vol = reconstruct(projs,rots,n_r,n_theta,max_shift,shift_step)

if ~exist('shift_step','var')
    shift_step = 0.5;
end

if ~exist('max_shift','var')
    max_shift = ceil(size(projs,1)*0.15); % max_shift is 15% of the image size
end

if ~exist('n_theta','var')
    n_theta = 360;
end

if ~exist('n_r','var')
    n_r = ceil(size(projs,1)*0.5);
end

[pf,~] = cryo_pft(projs,n_r,n_theta,'single');  % take Fourier transform of projections
pfC2 = cat(3,pf,pf); % Replicate the images

g = diag([-1 -1 1]);

nImages = size(rots,3);

RsC2 = zeros(3,3,2*nImages);

for k=1:nImages
    rot = rots(:,:,k);
    RsC2(:,:,k)           = rot;
    RsC2(:,:,nImages+k)   = g*rot;
end

[dxC2,~] = cryo_estimate_shifts(pfC2,RsC2,max_shift,shift_step,10000,[],0);

n = size(projs,1);
projsC2 = cat(3,projs,projs);
[ v1, ~, ~ ,~, ~, ~] = recon3d_firm( projsC2,RsC2,-dxC2, 1e-6, 100, zeros(n,n,n));
ii1=norm(imag(v1(:)))/norm(v1(:));
log_message('Relative norm of imaginary components = %e\n',ii1);
vol=real(v1);