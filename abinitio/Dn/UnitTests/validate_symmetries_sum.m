run('/scratch/yoelsh/aspire/initpath.m');
symmetry_degree = 40; % n of Dn.
rotations_num = 20; % Number of sampled rotations.
rotations = rand_rots(rotations_num);

angles = linspace(0, 360, symmetry_degree + 1);
angles = angles(2:symmetry_degree);
x_rot = quaternion([180, 0, 0], 'eulerd','XYZ','frame');
cn_symmetries = quaternion([zeros(numel(angles), 2), angles'], 'eulerd', 'XYZ', 'frame')';

for rotation_ind = 1:rotations_num
    rotation_i = rotations(:, :, rotation_ind).';  % Picking R^{T}_{i}
    for other_rotataion_ind = rotation_ind:rotations_num
        rotation_j = rotations(:, :, other_rotataion_ind); % Picking R_{j}
        
        % Transposing R_{j} is required only for the rotatepoint function.
        shifted_rotation_j = rotatepoint(x_rot, rotation_j.').';
        
        % Initializing averages to their first terms (where the symmetry is
        % the identity).
        rotations_average = rotation_i * rotation_j;
        other_symmetries_average = rotation_i * shifted_rotation_j;
        
        % Averaging over all symmetries of Cn.
        for cn_symmetry = cn_symmetries
            rotations_average = rotations_average + rotation_i * rotatepoint(cn_symmetry, rotation_j.').';
            other_symmetries_average = other_symmetries_average + rotation_i * rotatepoint(cn_symmetry, shifted_rotation_j.').';
        end
        rotation_average = rotations_average / symmetry_degree;
        other_symmetries_average = other_symmetries_average / symmetry_degree;
        
        % Asserting the Cn rotations average is indeed the outer-product of
        % the third rows of R_{i} and R_{j}.
        assert(all(abs(rotation_average - rotation_i(:, 3) * rotation_j(3, :)) < 1e-12, 'all'));
        
        % Asserting the other symmetries average is equal to the Cn
        % rotations average in magnitude, but different in sign.
        assert(all(abs(rotation_average + other_symmetries_average) < 1e-12, 'all'));
        
        % Assert the rotation average are all of rank 1.
        assert(rank(rotations_average) == 1);
    end
end
log_message('Validation Succeeded!');

