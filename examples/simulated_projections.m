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
rot_matrices = rand_rots(n);
% find inverse rotation matrices
inv_rot_matrices = permute(rot_matrices, [2 1 3]);

%% Generate 2D projections.
% Caclulate the 2D proejctions of the 3D volume in the directions
% determined by the randomly generated rotations matrices.

projections=cryo_project(ph,rot_matrices, 65,'single'); % generate phantom projecitons

% The following step is ESSENTIAL before calling FIRM reconstruction on the
% projections.
projections2=permute(projections,[2 1 3]);   % transpose each image


%% Reconstruct.
% Reconstruct the 3D volume from the 2D projections using FIRM.

[vol_recon, ~, ~ ,err, iter, flag] = recon3d_firm_parallel(projections2,inv_rot_matrices,[],1e-6,100,[]);
fprintf('Reconstruction error is %f\n', norm(ph(:)-vol_recon(:))/norm(ph(:))); 

% The printed error should about 0.004
