function [outputArg1,outputArg2] = estimate_simulation_error(rotations,...
    original_ests, estimated_rotations, symmetry_degree, theta_est,...
    order_2_gen)
%ESTIMATE_SIMULATION_ERROR Summary of this function goes here
%   Detailed explanation goes here
rotations_num = size(rotations, 3);

disp('Expected angles:');
print_relative_errors(rotations, original_ests, rotations_num, ...
    symmetry_degree);
disp('Output angles:');
disp(theta_est);
%         are_relative_rotations_correct(estimated_rotations(:,:,3), ...
%             rotations(:,:,3), estimated_rotations(:,:,4), rotations(:,:,4), ...
%             symmetry_degree);
% estimated_rotations_2 = estimate_inplane_rotations(...
%    pf, third_rows, symmetry_degree);

print_relative_errors(rotations, estimated_rotations, ...
    rotations_num, symmetry_degree);
dn_error = calc_estimation_error(rotations, estimated_rotations, ...
    symmetry_degree, order_2_gen);
fprintf('Total estimation error D%d w.r.t R (radians): %f\n', ...
    symmetry_degree, dn_error);

dof_O = quaternion([0, 0, 180 / symmetry_degree], 'eulerd', 'XYZ', 'frame')';
dof_O = quat2rotm(dof_O);
estimated_rotations2 = multiprod(dof_O, estimated_rotations, [1 2], [1 2]);

print_relative_errors(rotations, estimated_rotations2, rotations_num,...
    symmetry_degree);
dn_error2 = calc_estimation_error(rotations, estimated_rotations2, ...
    symmetry_degree, order_2_gen);
fprintf('Total estimation error D%d w.r.t OR (radians): %f\n', ...
    symmetry_degree, dn_error2);

dof_J = diag([1,1,-1]);
estimated_rotations3 = multiprod(dof_J, multiprod(estimated_rotations, ...
    dof_J, [1 2], [1 2]), [1 2], [1 2]);
print_relative_errors(rotations, estimated_rotations2, rotations_num,...
    symmetry_degree);
dn_error3 = calc_estimation_error(rotations, estimated_rotations3, ...
    symmetry_degree, order_2_gen);
fprintf('Total estimation error D%d w.r.t JRJ (radians): %f\n', ...
    symmetry_degree, dn_error3);

estimated_rotations4 = multiprod(dof_J, multiprod(estimated_rotations2, ...
    dof_J, [1 2], [1 2]), [1 2], [1 2]);
print_relative_errors(rotations, estimated_rotations2, rotations_num,...
    symmetry_degree);
dn_error4 = calc_estimation_error(rotations, estimated_rotations4, ...
    symmetry_degree, order_2_gen);
fprintf('Total estimation error D%d w.r.t OJRJ (radians): %f\n', ...
    symmetry_degree, dn_error4);

minimal_dn_error = min([dn_error, dn_error2, dn_error3, dn_error4]);
fprintf('Minimal estimation error D%d (radians): %f\n', ...
    symmetry_degree, minimal_dn_error);
end

function print_relative_errors(rotations, estimated_rotations,...
    rotations_num, n_symm)
    relative_rotations = multiprod(rotations, ...
        permute(estimated_rotations, [2 1 3]), [1 2], [1 2]);
    for i = 1:rotations_num
           fprintf('Relative Error %d (degrees):' , i);
           euler_angles = rad2deg(rotm2eul(relative_rotations(:, :, i),...
               'XYZ'));
           euler_angles(3) = mod(euler_angles(3), 360 / n_symm);
           disp(euler_angles);
    end
end

function rms_error = calc_estimation_error(rotations, ...
    estimated_rotations, symmetry_degree, order_2_gen)
rotations_num = size(rotations, 3);

angles = linspace(0, 360, symmetry_degree + 1);
angles = angles(1:symmetry_degree);
cn_symmetries = quat2rotm(quaternion([zeros(numel(angles), 2), angles'], ...
    'eulerd', 'XYZ', 'frame'));
gx_symmetries = multiprod(order_2_gen, cn_symmetries, [1 2], [1 2]);
dn_elements = cat(3, cn_symmetries, gx_symmetries);
dn_elements_count = 2 * symmetry_degree;
group_element_errors = NaN(dn_elements_count, rotations_num);

for k=1:rotations_num
    rot = rotations(:, :, k);
    estimation = estimated_rotations(:, :, k);
    
    for i = 1:dn_elements_count
        gn_power = dn_elements(:, :, i);
        moved_estimation = gn_power * estimation;
        
        group_element_errors(i, k) = geodesic_distance(rot, ...
            moved_estimation) ^ 2;
    end
end

[rotation_errors, I] = min(group_element_errors);
gn_min_elements = dn_elements(:, :, I);
registered_estimations = multiprod(gn_min_elements, ...
    estimated_rotations, [1 2], [1 2]);
%     for i = 1:rotations_num
%         fprintf('Ground Truth no. %d\n:' , i);
%         disp(rotations(:, :, i));
%         fprintf('Estimation no. %d\n:' , i);
%         disp(registered_estimations(:, :, i));
%     end
fprintf('Geodesic errors (degrees):');
disp(rad2deg(rotation_errors));
rms_error = sqrt(mean(rotation_errors));
end

