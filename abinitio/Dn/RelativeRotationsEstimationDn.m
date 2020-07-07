function RelativeRotationsEstimationDn(rotations, symmetry_degree)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

rotations_num = size(rotations, 3);  % N - number of rotations to estimate.

for rotation_ind = 1:rotations_num
    rotation_i = rotations(:, :, rotation_ind).';  % Picking R^{T}_{i}
    for other_rotataion_ind = rotation_ind:rotations_num
        rotation_j = rotations(:, :, other_rotataion_ind); % Picking R_{j}
        
        relative_rotation_cn = rotation_i(:, 3) * rotation_j(3, :);
    end
end
end
