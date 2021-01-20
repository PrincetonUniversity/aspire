run('/scratch/yoelsh/aspire/initpath.m');
addpath('../Utils', '..');
rotations_num = 20;
symmetry_degree = 3;  % n of Dn
seed = 1995;
rng(seed);
snr = 0;
noise_type = 'gaussian';
rotations = 0;
pf = 0;

if symmetry_degree == 4
    [rotations, pf] = generate_D4_molecule_data(rotations_num, snr, ...
        noise_type, seed);
elseif symmetry_degree == 3
    [rotations, pf] = generate_D3_phantom_data(rotations_num, snr, ...
        noise_type, seed);
else
    error("Symmetry degree n=%d is not supported", symmetry_degree);
end

third_rows = squeeze(rotations(3, :, :));
random_sign = 2 * randi([0 1], [1 rotations_num]) - 1; % Random +-1.
random_sign = repmat(random_sign, [3, 1]);
signed_third_rows = third_rows .* random_sign;

tic
[estimated_rotations, theta_est, original_ests] = investigate_inplane_rotation(pf, ...
    signed_third_rows, symmetry_degree);
duration = toc;

fprintf('Duration %f\n', duration);
disp('Expected angles:');
print_relative_errors(rotations, original_ests, rotations_num, ...
    symmetry_degree);
disp('Output angles:');
disp(theta_est);
% estimated_rotations_2 = estimate_inplane_rotations(...
%    pf, third_rows, symmetry_degree);

print_relative_errors(rotations, estimated_rotations, rotations_num,...
    symmetry_degree);
dn_error = calc_estimation_error(rotations, estimated_rotations, ...
    symmetry_degree);
fprintf('Total estimation error Dn first alternative: %f\n', dn_error);

dof = quaternion([0, 0, 180 / symmetry_degree], 'eulerd', 'XYZ', 'frame')';
dof = quat2rotm(dof);
estimated_rotations = multiprod(dof, estimated_rotations);

print_relative_errors(rotations, estimated_rotations, rotations_num,...
    symmetry_degree);
dn_error = calc_estimation_error(rotations, estimated_rotations, ...
    symmetry_degree);
fprintf('Total estimation error Dn second alternative: %f\n', dn_error);
    

% print_relative_rotations(rotations, estimated_rotations, rotations_num);
% cn_error = calc_estimation_error(rotations, estimated_rotations_2, ...
%    symmetry_degree);
% fprintf('Total estimation error Cn: %f\n', cn_error);


function print_relative_errors(rotations, estimated_rotations,...
    rotations_num, n_symm)
    relative_rotations = multiprod(rotations, ...
        permute(estimated_rotations, [2 1 3]), [1 2], [1 2]);
    for i = 1:rotations_num
          fprintf('Relative Error %d:' , i);
          euler_angles = rad2deg(rotm2eul(relative_rotations(:, :, i),...
              'XYZ'));
          euler_angles(3) = mod(euler_angles(3), 360 / n_symm);
          disp(euler_angles);
    end
end

function geo_dist = geodesic_distance(rot1, rot2)
    relative_rot = rot1 * rot2.';
    quat = rotm2quat(relative_rot);
    geo_dist = sqrt(2) * atan2(norm(quat(2:end)), quat(1));
end

function rms_error = calc_estimation_error(rotations, ...
    estimated_rotations, symmetry_degree)
    rotations_num = size(rotations, 3);
    
    angles = linspace(0, 360, symmetry_degree + 1);
    angles = angles(1:symmetry_degree);
    cn_symmetries = quaternion([zeros(numel(angles), 2), angles'], ...
        'eulerd', 'XYZ', 'frame')';
    gx_symmetries = quaternion([180*ones(numel(angles), 1), ...
        zeros(numel(angles), 1), angles'], 'eulerd', 'XYZ', 'frame')';
    dn_elements = cat(2, cn_symmetries, gx_symmetries);
    dn_elements_count = 2 * symmetry_degree;
    group_element_errors = NaN(dn_elements_count, rotations_num);
    
    for k=1:rotations_num
        rot = rotations(:, :, k);
        estimation = estimated_rotations(:, :, k);
        
        for i = 1:dn_elements_count
            gn_power = dn_elements(1, i);
            true_rot = rotatepoint(gn_power, rot.').';
            
            group_element_errors(i, k) = geodesic_distance(true_rot, ...
                estimation) ^ 2;       
        end
    end
    
    rotation_errors = min(group_element_errors);
    fprintf('Geodesic errors:');
    disp(rad2deg(rotation_errors));
    rms_error = sqrt(mean(rotation_errors));
end

function rr = get_relative_rotations(rotations, n_symm)
    rotations_num = size(rotations, 3);
    rr = cell(rotations_num, rotations_num);
    angles = linspace(0, 360, n_symm + 1);
    angles = angles(1:n_symm);
    cn_symmetries = quat2rotm(quaternion([zeros(numel(angles), 2), angles'], ...
        'eulerd', 'XYZ', 'frame')');
    gx_symmetries = quat2rotm(quaternion([180*ones(numel(angles), 1), ...
        zeros(numel(angles), 1), angles'], 'eulerd', 'XYZ', 'frame')');
    dn_elements = cat(3, cn_symmetries, gx_symmetries);
    
    for i = 1:rotations_num
        for j = 1:rotations_num
            rr{i, j} = multiprod(multiprod(rotations(:, :, i).', dn_elements), ...
                rotations(:, :, j));
        end
    end
end
