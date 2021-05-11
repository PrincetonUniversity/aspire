function [estimated_rotations, theta_est, original_ests] = runDn(...
    pf_images, symmetry_degree, gx, varargin)
% RUNDN computes a reconstruction of a D2 volume from a given set of
%   projection images.
%
%   Input:
%       pf_images - Polar Fourier-Transformed projection images. The i-th
%                   image is given by pf_images(:,:,i).
%       symmetry_degree - an integer >= 3.
%       pf_images - Polar Fourier-Transformed projection images. The i-th
%                   image is given by pf_images(:,:,i).
%       gx - Either the matrix diag([1,-1,-1]) or diag([-1,1,-1]).
%
%   Output:
%       TBD
%       
%   Written by Elad Eatah April 2021.

rotations_num = size(pf_images, 3);
group_elements = cryo_Dn_group_elements(symmetry_degree, gx);
rotations = 0;
if nargin == 4
    rotations = varargin{1};
end

%% Estimate relative-rotations between all desired rotations.
relative_rotations = relative_rotation_estimation_Dn(symmetry_degree, ...
    pf_images, group_elements);%, rotations);

%% Handedness-Synchronization stage
% Estimates outer-products of third-rows, which either all have spurios J
% or none of them have such J.
vijs = HandednessSynchronizationDn(relative_rotations, rotations_num, ...
    symmetry_degree);

%% Estimate third-rows, up to sign, of the desired rotations.
[signed_third_rows, ~, ~] = syncSignsDn(vijs, rotations_num);
signed_third_rows = signed_third_rows';

%% Inplane rotation estimation
[estimated_rotations, theta_est, original_ests] = ...
    InPlaneRotationEstimation(pf_images, signed_third_rows, ...
    symmetry_degree, gx);
end
