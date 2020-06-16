function RelativeRotationsEstimationDn(rotations, symmetry_degree)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

for rotation_ind = 1:rotations_num
    rotation_i = rotations(:, :, rotation_ind).';  % Picking R^{T}_{i}
    for other_rotataion_ind = rotation_ind:rotations_num
        rotation_j = rotations(:, :, other_rotataion_ind); % Picking R_{j}
        
        relative_rotation_cn = rotation_i(:, 3) * rotation_j(3, :);
        relative_rotation_compliment = -relative_rotation_cn;
    end
end
end
