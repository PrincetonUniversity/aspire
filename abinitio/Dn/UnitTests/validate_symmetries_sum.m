run('/scratch/yoelsh/aspire/initpath.m');
symmetry_degree = 16; % n of Dn.
rotations_num = 20; % Number of sampled rotations.
rotations = rand_rots(rotations_num);

angles = linspace(0, 360, symmetry_degree + 1);
angles = angles(2:symmetry_degree);
x_rot = quaternion([180, 0, 0], 'eulerd','XYZ','frame');
y_rot = quaternion([0, 180, 0], 'eulerd','XYZ','frame');
cn_symmetries = quaternion([zeros(numel(angles), 2), angles'], 'eulerd', 'XYZ', 'frame')';

for rotation_ind = 1:rotations_num
    % Picking R_{i} and its transpose.
    rotation_i = rotations(:, :, rotation_ind);
    rotation_i_T = rotation_i.';
    % Transposing here is required only for the rotatepoint function.
    shifted_rotation_i = rotatepoint(x_rot, rotation_i_T);
    
    shifted_rotation_i_rot_y = rotatepoint(y_rot, rotation_i_T);
    
    for other_rotataion_ind = rotation_ind:rotations_num
        rotation_j = rotations(:, :, other_rotataion_ind); % Picking R_{j}
        
        % Initializing averages to their first terms (where the symmetry is
        % the identity).
        rotations_average = rotation_i_T * rotation_j;
        other_symmetries_average = shifted_rotation_i * rotation_j;
        
        % Assert claim 0.4 for the specific case of s=0
        mean_relative_rotation = (rotations_average + other_symmetries_average) / 2;
        assert(all(abs(mean_relative_rotation - rotation_i(1, :).' * rotation_j(1, :)) < 1e-12, 'all'));
        weighted_rotations_mean = 2 * mean_relative_rotation;
        
        % Assert second rows outer product identity for the specific case
        % of s = 0
        weighted_rotations_mean_rot_y = (rotations_average + shifted_rotation_i_rot_y * rotation_j) / 2;
        assert(all(abs(weighted_rotations_mean_rot_y - rotation_i(2, :).' * rotation_j(2, :)) < 1e-12, 'all'));
        weighted_rotations_mean_rot_y = 2 * weighted_rotations_mean_rot_y;
        n_s = symmetry_degree - 1;
        
        % Averaging over all symmetries of Cn.
        for cn_symmetry = cn_symmetries
            shifted_rotation_j = rotatepoint(cn_symmetry, rotation_j.').';
            
            current_rotation = rotation_i_T * shifted_rotation_j;
            rotations_average = rotations_average + current_rotation;            
            
            current_gx_rotation = shifted_rotation_i * shifted_rotation_j;
            other_symmetries_average = other_symmetries_average + current_gx_rotation;
            
            % Assert claim 0.4 for any s >= 1
            mean_relative_rotation = (current_rotation + current_gx_rotation) / 2;
            rot_line = quat2rotm(cn_symmetry);
            first_rot_line = rot_line(1, :);
            assert(all(abs(mean_relative_rotation - rotation_i(1, :).' * first_rot_line * rotation_j) < 1e-12, 'all'));            
            assert(all(abs(mean_relative_rotation - cos(2 * pi * n_s / symmetry_degree) ...
                * rotation_i(1, :).' * rotation_j(1, :) + sin(2 * pi * n_s / symmetry_degree) ...
                * rotation_i(1, :).' * rotation_j(2, :)) < 1e-12, 'all'));
            
            current_gy_rotation = shifted_rotation_i_rot_y * shifted_rotation_j;
            mean_relative_rot_y = (current_rotation + current_gy_rotation) / 2;
            second_rot_line = rot_line(2, :);
            assert(all(abs(mean_relative_rot_y - rotation_i(2, :).' * second_rot_line * rotation_j) < 1e-12, 'all'));
            assert(all(abs(mean_relative_rot_y - sin(2 * pi * n_s / symmetry_degree) ...
                * rotation_i(2, :).' * rotation_j(1, :) - cos(2 * pi * n_s / symmetry_degree) ...
                * rotation_i(2, :).' * rotation_j(2, :)) < 1e-12, 'all'));
            
            
            weighted_rotations_mean = weighted_rotations_mean + ...
                cos(2 * pi * n_s / symmetry_degree) * (current_rotation + current_gx_rotation);
            weighted_rotations_mean_rot_y = weighted_rotations_mean_rot_y +...
                cos(2 * pi * n_s / symmetry_degree) * (current_rotation + current_gy_rotation);
            n_s = n_s - 1;
        end
        rotation_average = rotations_average / symmetry_degree;
        other_symmetries_average = other_symmetries_average / symmetry_degree;
        
        weighted_rotations_mean = weighted_rotations_mean / symmetry_degree;
        assert(all(abs(weighted_rotations_mean - rotation_i(1, :).' * rotation_j(1, :)) < 1e-12, 'all'));
        
        weighted_rotations_mean_rot_y = weighted_rotations_mean_rot_y / symmetry_degree;
        assert(all(abs(weighted_rotations_mean_rot_y - rotation_i(2, :).' * rotation_j(2, :)) < 1e-12, 'all'));
        
        % Asserting the Cn rotations average is indeed the outer-product of
        % the third rows of R_{i} and R_{j}.
        assert(all(abs(rotation_average - rotation_i(3, :).' * rotation_j(3, :)) < 1e-12, 'all'));
        
        % Asserting the other symmetries average is equal to the Cn
        % rotations average in magnitude, but different in sign.
        assert(all(abs(rotation_average + other_symmetries_average) < 1e-12, 'all'));
        
        % Assert the rotation average are all of rank 1.
        assert(rank(rotations_average) == 1);
    end
end
log_message('Validation Succeeded!');
clear;

