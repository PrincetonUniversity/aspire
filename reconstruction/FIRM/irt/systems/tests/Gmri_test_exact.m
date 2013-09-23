% Gmri_test_exact.m
% test "exact" option of Gmri object

if 1 % test integer k-space samples and unit fov (for KK)
	% todo: add warning to nufft_init if DFT sampling with small nbhd
	% and KB interpolator
	kspace = [-3:2:1 0.5]'; kspace = [kspace ones(size(kspace))];
	N = [16 8];
	Ge = Gmri(kspace, true(N), 'exact', true, 'n_shift', N/2);
	Ge = Ge(:,:);
	G0 = Gmri(kspace, true(N));
	G0 = G0(:,:);
	minmax(Ge-G0)
	G1 = Gmri(kspace, true(N), 'nufft', ...
		{N, [1 1], 1*N, N/2, 'uniform'});
%		{N, [2 2], 1*N, N/2, 'table', 2^10, 'minmax:kb'}); % big error!
	G1 = G1(:,:);
	minmax(Ge-G1)
%	G1 = Gmri(kspace, true(N), 'nufft', {
end

ig = image_geom('nx', 64, 'ny', 60, 'dx', 1);
M = 100;
ti = linspace(0, 30e-3, M); % 30 msec readout
zmap = 20 * ig.ones + 2i * pi * 10;
kspace = zeros(M,2); % FID only
%n_shift = [0 0];
G = Gmri(kspace, ig.mask, 'exact', true, ...
	'ti', ti, 'zmap', zmap, 'n_shift', [0 0]);

tmp = feval(G.arg.new_image_basis, G, {'dirac'});

x = ig.ones;
x = x / sum(x(:));
yb = G * x(ig.mask);
if im
	plot(ti, real(yb), '-', ti, imag(yb), '--', ti, abs(yb), 'r-')
%	plot(ti, angle(yb), '-')
end
pr exp(-max(ti)*max(real(zmap(:))))

xb = G' * yb;
xb = ig.embed(xb);

if 1 % 1d test
	fov = 7;
	N = 128;
	dx = fov / N;
	n_shift = 17;
	r = ([0:N-1] - n_shift) * dx;
	k = [-N/2:N/2-0]' / fov;
	ti = linspace(0, 5e-3, length(k))'; % 5 ms
	zmap = 10 * [1:N]'/N + 15i; % 1/s and rad/s
	G0 = exp(-2i * pi * k * r(:)') .* exp(-ti * zmap.');
	G1 = Gmri(k, true(N,1), ...
		'ti', ti, 'zmap', zmap, 'L', 10, ...
		'fov', fov, 'exact', true, 'n_shift', n_shift);
	G1 = G1(:,:);
%	im_toggle(real(G0), real(G1))
%	equivs(real(G0), real(G1))
	max_percent_diff(G0, G1)
end
