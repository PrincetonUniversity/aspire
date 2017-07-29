%% Using simulated projections
%
% Usage example of functions for generating simulated projections.
%
% The script generates a volume sampled from a sum of 3D Gaussians,
% computes 2D projections of the sampled volume, and reconstructs back the
% volume using FIRM.
%
% Gene Katsevich and Yoel Shkolnisky, March 2014.

%% Generate volume.
% Generate a volume of size 65x65x65 sampled from a non-symmetric sum of
% Guassians   

ph=cryo_gaussian_phantom_3d('C1_params',65,1); % generate a 3D phantom

%% Generate random rotations.
% Create 100 random uniformly sampled rotations (from the uniform
% distribution on SO(3)).

n=100;                              % use 100 images
inv_rot_matrices = zeros(3, 3, n);
q = qrand(n);                       % generate rotations as quaternions
% find inverse rotation matrices
for k = 1:n
    quat = q(:, k);
    q(:, k) = quat/norm(quat);
    rot = q_to_rot(q(:, k));
    inv_rot_matrices(:, :, k) = rot';
end

%% Generate 2D projections.
% Caclulate the 2D proejctions of the 3D volume in the directions
% determined by the randomly generated rotations matrices.

projections=cryo_project(ph,q, 65,'single'); % generate phantom projecitons

% The following step is ESSENTIAL before calling FIRM reconstruction on the
% projections.
projections2=permute(projections,[2 1 3]);   % transpose each image


%% Reconstruct.
% Reconstruct the 3D volume from the 2D projections using FIRM.

params = struct();
params.rot_matrices = inv_rot_matrices;
params.ctf = ones(size(projections2, 1)*ones(1, 2));
params.ctf_idx = ones(n, 1);
params.shifts = zeros(2, n);
params.ampl = ones(n, 1);

basis = dirac_basis(size(projections2, 1)*ones(1, 3));

vol_recon = cryo_estimate_mean(projections2, params, basis);
fprintf('Reconstruction error is %f\n', norm(ph(:)-vol_recon(:))/norm(ph(:))); 

% The printed error should about 0.004
