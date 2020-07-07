run('/scratch/yoelsh/aspire/initpath.m');
symmetry_degree = 24; % n of Dn.
rotations_num = 100; % Number of sampled rotations.
tol = 1e-12;
rotations = rand_rots(rotations_num);
third_rows = squeeze(rotations(3, :, :)).';
relative_rotations = zeros(3, 3, nchoosek(rotations_num, 2));

angles = linspace(0, 360, symmetry_degree + 1);
angles = angles(2:symmetry_degree);
cn_symmetries = quaternion([zeros(numel(angles), 2), angles'], 'eulerd', 'XYZ', 'frame')';

for rotation_ind = 1:rotations_num
    rotation_i_transpose = rotations(:, :, rotation_ind).';  % Picking R^{T}_{i}
    for other_rotation_ind = rotation_ind+1:rotations_num
        rotation_j = rotations(:, :, other_rotation_ind); % Picking R_{j}
        
        % Initializing averages to their first terms (where the symmetry is
        % the identity).
        rotations_average = rotation_i_transpose * rotation_j;
        
        % Averaging over all symmetries of Cn.
        for cn_symmetry = cn_symmetries
            rotations_average = rotations_average + rotation_i_transpose * rotatepoint(cn_symmetry, rotation_j.').';
        end
        current_relative_rotation = rotations_average / symmetry_degree;
        
        % Random switching between mean over Cn symmetries and Dn-Cn
        % symmetries.
        if rand <= 0.5
            current_relative_rotation = -current_relative_rotation;
        end
        
        current_pair_index = uppertri_ijtoind(rotation_ind, other_rotation_ind, rotations_num);
        relative_rotations(:, :, current_pair_index) = current_relative_rotation;
    end
end

addpath('../', '../D2/Synchronization');
% Adapting variables for Eitan's script for D2.
cVec = repmat([1, 2, 3], 1, nchoosek(rotations_num, 2));
relative_rotations_d2 = repmat(relative_rotations, 1, 1, 1, 3);
[U1_d2, a_d2, b_d2] = syncSigns(relative_rotations_d2, cVec, rotations_num); % This call is OK.
U1_d2 = squeeze(U1_d2(1, :, :))';

% The adapted version for sign synchronization.
[U1_dn, a_dn, b_dn] = syncSignsDn(relative_rotations, rotations_num);

is_d2_sign_symmetry_ok = abs(U1_d2 - third_rows) < tol | abs(U1_d2 + third_rows) < tol;
is_dn_sign_symmetry_ok = abs(U1_dn - third_rows) < tol | abs(U1_dn + third_rows) < tol;

assert(all(is_d2_sign_symmetry_ok, 'all'));
assert(all(is_dn_sign_symmetry_ok, 'all'));
log_message('Sign Synchronization for Dn is successfull!');