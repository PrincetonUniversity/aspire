function vol = reconstruct_ml_c12(projs,rots,n_r,n_theta,max_shift,shift_step)

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
pfC12 = cat(3,pf,pf,pf,pf,pf,pf,pf,pf,pf,pf,pf,pf); % Replicate the images

g = [cosd(360/12) -sind(360/12) 0; 
	 sind(360/12)  cosd(360/12) 0; 
	 0 				 0  1]; % a rotation of 90 degrees about the z-axis

nImages = size(rots,3);

RsC12 = zeros(3,3,12*nImages);

for k=1:nImages
    rot = rots(:,:,k);
    for s=0:11
        RsC12(:,:,s*nImages+k) = g^s*rot; 
    end
end

[dxC12,~] = cryo_estimate_shifts(pfC12,RsC12,ceil(2*sqrt(2)*max_shift),shift_step,10000,[],0);

n = size(projs,1);
projsC12 = cat(3,projs,projs,projs,projs,projs,projs,projs,projs,projs,projs,projs,projs);
[ v1, ~, ~ ,~, ~, ~] = recon3d_firm( projsC12,RsC12,-dxC12, 1e-6, 100, zeros(n,n,n));
ii1=norm(imag(v1(:)))/norm(v1(:));
log_message('Relative norm of imaginary components = %e\n',ii1);
vol=real(v1);