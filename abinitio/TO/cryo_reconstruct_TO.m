function [vol,rotations, est_shifts]= cryo_reconstruct_TO(symmetry,projs,rots,n_r,n_theta,max_shift,shift_step)

% cryo_reconstruct_TOreconsturction of a T or O symmetric molecule
%
% Parameters
%   symmetry          symmetry type: 'T' or 'O'.
%   projs             array of projection images. 
%   rots              array of estimated rotation matrices for each projection image.
%   n_theta           Angular resolution for common lines detection.
%                     Default 360.
%   n_r               Radial resolution for common line detection as a
%                     percentage of image size. Default is half the width of the images.
%   max_shift         Maximal 1d shift (in pixels) to search between
%                     common-lines. Default is 15% of image width of the images.
%   shift_step        Resolution of shift estimation in pixels. Note
%                     that shift_step can be any positive real number. Default:0.5.
%
%   Written by Adi Shasha January 2021. 

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

gR = cryo_TO_group_elements(symmetry);
n_gR = size(gR,3);

[pf,~] = cryo_pft(projs,n_r,n_theta,'single');  % take Fourier transform of projections
pf_TO = [];
for i=1:n_gR
    pf_TO = cat(3,pf_TO,pf);
end

nImages = size(rots,3);

Rs_TO = zeros(3,3,nImages*n_gR);

for k=1:nImages
    rot = rots(:,:,k);
    for s=0:n_gR-1
        Rs_TO(:,:,s*nImages+k) = gR(:,:,s+1)*rot; 
    end
end


log_message('Estimating shifts');
[dx_TO,~] = cryo_estimate_shifts(pf_TO,Rs_TO,ceil(2*sqrt(2)*max_shift),shift_step,10000,[],0);
n = size(projs,1);
projs_TO = [];
for i=1:n_gR
    projs_TO = cat(3,projs_TO,projs);
end

log_message('Reconstructing');
%[ v1, ~, ~ ,~, ~, ~] = recon3d_firm( projsCn,RsCn,-dxCn, 1e-6, 100, zeros(n,n,n));
params = struct();
params.rot_matrices = Rs_TO;
params.ctf = ones(n*ones(1, 2));
params.ctf_idx = ones(size(projs_TO, 3), 1);
params.shifts = dx_TO.';
params.ampl = ones(size(projs_TO, 3), 1);

mean_est_opt.max_iter = 100;
mean_est_opt.rel_tolerance = 1.0e-6;
mean_est_opt.verbose = false;
mean_est_opt.precision = 'single';

basis = dirac_basis(size(projs_TO, 1)*ones(1, 3));

% v1 = cryo_estimate_mean(single(projs_TOI), params, basis);
v1 = cryo_estimate_mean(single(projs_TO), params, basis, mean_est_opt);

% ii1=norm(imag(v1(:)))/norm(v1(:));
% log_message('Relative norm of imaginary components = %e\n',ii1);
vol=real(v1);
% vol=v1;

rotations=Rs_TO;
est_shifts=dx_TO;
