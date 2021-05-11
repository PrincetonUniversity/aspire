%%  RUN_DN_MOLECULE_SIMULATION 
%   Performs some experiments of the entire algorithm for Dn-structures.
%   Each experiment is applied on a specific value of symmetry-degree
%   (n>2) and a specific value of SNR.

%% Configuration.
run('/scratch/yoelsh/aspire/initpath.m');
addpath('../Utils', '..');
gx = diag([1, -1, -1]);
gy = diag([-1, 1, -1]);
rotations_num = 10;  % Number of orientations to estimate (N).
symmetry_degrees = [3];  % n of Dn, a list of values for experiments.
order_2_gen = gx;  % X-axis is the 'first' symmetry axis in XY-plane.

use_yoels = false;
seed = 1995;  % Setting seed for reproducibility.
rng(seed);
snrs = [0];  % List of tested SNR values (non-negative or Inf).
noise_type = 'gaussian';  % Applied noise-type.
rotations = 0;
pf = 0;

%% Performing all experiments.
for symmetry_degree = symmetry_degrees
    for snr = snrs
        % Directing the log of this experiment to a text file no disk.
        log_file_name = sprintf('log D%d SNR %g.txt', symmetry_degree, snr);
        log_file_name = fullfile('Simulations', 'simulation_logs', ...
            log_file_name);
        diary (log_file_name)
        
        % Generate the simulation's input - Noisy projections of a 
        % Dn-structure and the orientations the algorithm seeks 
        % (ground-truth).
        if symmetry_degree > 2
            [rotations, pf] = generate_Dn_molecule_data(symmetry_degree, ...
                rotations_num, order_2_gen, snr, noise_type, seed, use_yoels);
        else
            error("Symmetry degree n=%d is not supported", symmetry_degree);
        end
        
        % Apply the entire algorithm for Dn.
        [estimated_rotations,~,original_ests] = runDn(pf, ...
            symmetry_degree, order_2_gen, rotations);
        
        % Display the errors of this simulation.
        estimate_simulation_error(rotations, original_ests, ...
            estimated_rotations, symmetry_degree, theta_est, order_2_gen);
        diary off
    end
end

function rr = get_relative_rotations(rotations, n_symm, order_2_gen)
    rotations_num = size(rotations, 3);
    rr = cell(rotations_num, rotations_num);
    angles = linspace(0, 360, n_symm + 1);
    angles = angles(1:n_symm);
    cn_symmetries = quat2rotm(quaternion([zeros(numel(angles), 2), angles'], ...
        'eulerd', 'XYZ', 'frame')');
    gx_symmetries = multiprod(order_2_gen, cn_symmetries);
    dn_elements = cat(3, cn_symmetries, gx_symmetries);
    
    for i = 1:rotations_num
        for j = 1:rotations_num
            rr{i, j} = multiprod(multiprod(rotations(:, :, i).', dn_elements), ...
                rotations(:, :, j));
        end
    end
end

function mean_rr_err = are_relative_rotations_correct(Ri_estimate, ...
    Ri_correct, Rj_estimate, Rj_correct, n_symm, order_2_gen)
    RR = NaN(3,3,2);
    RR(:,:,1) = Ri_estimate;
    RR(:,:,2) = Rj_estimate;
    rr = get_relative_rotations(RR, n_symm, order_2_gen);
    rr = rr{1, 2};
    RR(:,:,1) = Ri_correct;
    RR(:,:,2) = Rj_correct;
    correct_rr = get_relative_rotations(RR, n_symm, order_2_gen);
    correct_rr = correct_rr{1,2};
    %gx = diag([1, -1, -1]);
    mean_rr_err = 0;
    
    for i=1:2*n_symm
        min_err = 1e+17;
        correct_relative_rotation = correct_rr(:,:,i);
        for j=1:2*n_symm
            current_err = rad2deg(geodesic_distance(rr(:, :, j), ...
                correct_relative_rotation));
            min_err = min([min_err, current_err]);                
        end
        disp("Relative Rotation Error:")
        disp(min_err);
        mean_rr_err = mean_rr_err + min_err;
    end
    mean_rr_err = mean_rr_err / (2 * n_symm);
end

function rr_error = get_rr_error(estimated_rr, ground_truth, ...
    symmetry_degree, order_2_gen)
total_err = 0;
pairs_count = size(estimated_rr, 1);
rotations_num = size(ground_truth, 3);
correct_rr = get_relative_rotations(ground_truth, symmetry_degree, ...
    order_2_gen);
p = [];
pi = [];
pj = [];

for i_ind = 1:rotations_num
    for j_ind = i_ind+1:rotations_num
        ii = uppertri_ijtoind_vec(i_ind, j_ind, rotations_num);
        current_estimated_rr = squeeze(estimated_rr(ii,:,:,:));
        %[i_ind, j_ind] = ind2sub([rotations_num, rotations_num], ii);
        rr = correct_rr{i_ind, j_ind};
        
        mean_err = 0;
        for k=1:2*symmetry_degree
            min_err = 1e+17;
            correct_relative_rotation = rr(:,:,k);
            
            for j=1:2*symmetry_degree
                estimation_of_rr = current_estimated_rr(:, :, j);
                current_err = rad2deg(geodesic_distance(estimation_of_rr, ...
                    correct_relative_rotation));
                min_err = min([min_err, current_err]);
            end
            mean_err = mean_err + min_err;
            if min_err > 0.1
                p = [p, ii];
                pi = [pi, i_ind];
                pj = [pj, j_ind];
            end
        end
        fprintf("Relative Rotation Error pair no. %d (degrees):", ii);
        mean_err = mean_err / (2 * symmetry_degree);
        disp(mean_err);
        
        total_err = total_err + mean_err;
    end
end
total_err = total_err / pairs_count;
disp("Mean Relative Rotation Error (degrees):");
disp(total_err);
disp("Problematic indices:");
disp(unique(p));
disp("Problematic indices ij :");
disp(unique([pi; pj]', 'rows')');

end