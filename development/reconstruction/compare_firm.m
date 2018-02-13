n = 65;
K = 1000;
SNR = 1000;
max_shift = 0.1*n;
shift_step = max_shift/10;
defocus = linspace(1000, 2500, 7);
[projections, ~, shifts, rotations] = cryo_gen_projections(n, K, SNR, ...
    max_shift, shift_step);

ctf = zeros([n*ones(1, 2) numel(defocus)]);
ctf_idx = mod(0:K-1, numel(defocus))+1;

for k = 1:numel(defocus)
    ctf(:,:,k) = cryo_CTF_Relion(n, 300, defocus(k), defocus(k), 0, 2.0, ...
        1.4, 0.1);

    projections(:,:,ctf_idx == k) = ...
        cryo_addctf(projections(:,:,ctf_idx == k), ctf(:,:,k));
end

tol = 1e-6;
max_iter = 1000;

inv_rotations = invert_rots(rotations);

fprintf('Reconstruction from clean projections\n');
tic;
v1 = recon3d_firm_ctf(projections, ctf, ctf_idx, inv_rotations, -shifts, ...
    tol, max_iter, []);
t1 = toc;

f = load(fullfile('projections', 'simulation', 'maps', 'cleanrib.mat'));
volref = f.volref;

fprintf('Reconstruction took %.2f seconds\n', t1);
fprintf('The relative error of reconstruction is %f.\n',...
    norm(v1(:)-volref(:))/norm(volref(:)));

params = struct();
params.rot_matrices = inv_rotations;
params.ctf = ctf;
params.ctf_idx = ctf_idx;
params.shifts = shifts.';
params.ampl = ones(size(projections, 3), 1);

mean_est_opt.max_iter = max_iter;
mean_est_opt.rel_tolerance = tol;
mean_est_opt.verbose = false;
mean_est_opt.precision = 'single';

basis = dirac_basis(size(projections, 1)*ones(1, 3));

tic;
v2 = cryo_estimate_mean(single(projections), params, basis, mean_est_opt);
t2 = toc;

fprintf('Reconstruction took %.2f seconds\n', t2);
fprintf('The relative error of reconstruction is %f.\n',...
    norm(v2(:)-volref(:))/norm(volref(:)));
