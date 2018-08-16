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
pfC4 = cat(3,pf,pf,pf,pf); % Replicate the images

g = [0 -1 0; 1 0 0; 0 0 1]; % a rotation of 90 degrees about the z-axis

nImages = size(rots,3);

RsC4 = zeros(3,3,4*nImages);

for k=1:nImages
    rot = rots(:,:,k);
    RsC4(:,:,k)           = rot;
    RsC4(:,:,nImages+k)   = g*rot;
    RsC4(:,:,2*nImages+k) = g*g*rot;
    RsC4(:,:,3*nImages+k) = g*g*g*rot;
end

[dxC4,~] = cryo_estimate_shifts(pfC4,RsC4,max_shift,shift_step,10000,[],0);

n = size(projs,1);
projsC4 = cat(3,projs,projs,projs,projs);
%[ v1, ~, ~ ,~, ~, ~] = recon3d_firm( projsC4,RsC4,-dxC4, 1e-6, 100, zeros(n,n,n));
params = struct();
params.rot_matrices = RsC4;
params.ctf = ones(n*ones(1, 2));
params.ctf_idx = ones(size(projsC4, 3), 1);
params.shifts = dxC4.';
params.ampl = ones(size(projsC4, 3), 1);
mean_est_opt.max_iter = 100;
mean_est_opt.rel_tolerance = 1.0e-6;
mean_est_opt.verbose = false;
mean_est_opt.precision = 'single';
basis = dirac_basis(size(projsC4, 1)*ones(1, 3));
v1 = cryo_estimate_mean(single(projsC4), params, basis, mean_est_opt);

ii1=norm(imag(v1(:)))/norm(v1(:));
log_message('Relative norm of imaginary components = %e\n',ii1);
vol=real(v1);