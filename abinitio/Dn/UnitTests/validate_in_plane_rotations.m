%%  VALIDATE_RELATIVE_ROTATION_ESTIMATION
%   Evaluates the performance of the relative-rotations estimation step of
%   the algorithm for Dn. When no noise is applied, and the candidates set
%   equals the ground-truth matrices, the relative-rotations should be
%   estimated correctly.

%% Configuration
run('/scratch/yoelsh/aspire/initpath.m');
addpath('../Utils', '..');
gx = diag([1, -1, -1]);
gy = diag([-1, 1, -1]);
rotations_num = 10;
symmetry_degrees = [3];  % n of Dn
order_2_gen = gx;  % X-axis is the 'first' symmetry axis in XY-plane.

use_yoels = false;
seed = 1995;
rng(seed);
snrs = [0];
noise_type = 'gaussian';
rotations = 0;
pf = 0;

%% Performing all experiments.
for symmetry_degree = symmetry_degrees
    for snr = snrs
        log_file_name = sprintf("inplane rotation estimation " + ...
            "log D%d SNR %g.txt", symmetry_degree, snr);
        log_file_name = fullfile('Simulations', 'simulation_logs', ...
            log_file_name);
        diary (log_file_name)
        
        if symmetry_degree > 2
            [rotations, pf] = generate_Dn_molecule_data(symmetry_degree, ...
                rotations_num, order_2_gen, snr, noise_type, seed, use_yoels);
        else
            error("Symmetry degree n=%d is not supported", symmetry_degree);
        end
        
        signed_third_rows = pick_ground_truth_third_rows(rotations);
        
        tic
        [estimated_rotations, theta_est, original_ests] = investigate_inplane_rotation(pf, ...
            signed_third_rows, symmetry_degree, gx);
        duration = toc;

        fprintf('Duration %f\n', duration);
        estimate_simulation_error(rotations, original_ests, ...
            estimated_rotations, symmetry_degree, theta_est, order_2_gen);
        diary off
    end
end

%% Function for picking third-rows, with random signs
function signed_third_rows = pick_ground_truth_third_rows(rotations)
    rotations_num = size(rotations, 3);
    third_rows = squeeze(rotations(3, :, :));
    random_sign = 2 * randi([0 1], [1 rotations_num]) - 1; % Random +-1.
    random_sign = repmat(random_sign, [3, 1]);
    signed_third_rows = third_rows .* random_sign;
end