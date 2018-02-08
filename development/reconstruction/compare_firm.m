n = 17;
K = 1000;
SNR = 1000;
[projections, ~, ~, rotations] = cryo_gen_projections(n, K, SNR);

tol = 1e-6;
max_iter = 1000;

inv_rotations = invert_rots(rotations);

fprintf('Reconstruction from clean centered projections\n');
tic;
v1 = recon3d_firm(projections, inv_rotations, [], tol, max_iter, []);
t1 = toc;

f = load(fullfile('projections', 'simulation', 'maps', 'cleanrib.mat'));
volref = f.volref;

fprintf('Reconstruction took %.2f seconds\n', t1);
%fprintf('The relative error of reconstruction is %f.\n',...
    %norm(v1(:)-volref(:))/norm(volref(:)));

params = struct();
params.rot_matrices = inv_rotations;
params.ctf = ones(size(projections, 1)*ones(1, 2));
params.ctf_idx = ones(size(projections, 3), 1);
params.shifts = zeros(2, size(projections, 3));
params.ampl = ones(size(projections, 3), 1);

mean_est_opt.max_iter = max_iter;
mean_est_opt.rel_tolerance = tol;
mean_est_opt.verbose = false;

basis = dirac_basis(size(projections, 1)*ones(1, 3));

tic;
v2 = cryo_estimate_mean(projections, params, basis, mean_est_opt);
t2 = toc;

fprintf('Reconstruction took %.2f seconds\n', t2);
%fprintf('The relative error of reconstruction is %f.\n',...
    %norm(v2(:)-volref(:))/norm(volref(:)));
