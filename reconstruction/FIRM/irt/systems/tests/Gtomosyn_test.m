% Gtomosyn_test.m
% Copyright 2010-07-27, Jeff Fessler, University of Michigan

% small system for basic test
if ~isvar('As'), printm 'setup small'
	f.down = 16;
	f.dx = 0.1; % 100 um
	f.dz = 1; % 1 mm
	igs = image_geom('nx', 2304, 'ny', 48, 'nz', 1920, ...
		'dx', f.dx, 'dy', f.dz, 'dz', f.dx, ...
		'down', f.down);

	cgs = ct_geom('fan', 'ns', 2304, 'nt', 1920, 'na', 1, ...
			'orbit', 0, ...
			'ds', f.dx, 'dt', f.dx, ...
			'dfs', inf, ... % flat
			'dsd', 650, ...
			'dod', abs(igs.ny*igs.dy/2), ... % compressed on detector
			'down', f.down);

	As = Gtomosyn(cgs, igs, 'chat', 0); % set to 1 to see geometry!

	tmp = As * igs.unitv;
	im(tmp), cbar

%	Fatrix_test_basic(As, igs.mask)
	tester_tomo2(As, igs.mask, 'halt', 0, 'nblock', 0) % todo: ordering problem

	test_adjoint(As, 'big', 1, 'tol', 2e-6)
prompt
end


% big system for timing test
if ~isvar('A'), printm 'setup big'
	f.down = 4;
	f.dx = 0.1; % 100 um
	f.dz = 1; % 1 mm
	ig = image_geom('nx', 2304, 'ny', 48, 'nz', 1920, ...
		'dx', f.dx, 'dy', f.dz, 'dz', f.dx, ...
		'down', f.down);

	cg = ct_geom('fan', 'ns', 2304, 'nt', 1920, 'na', 1, ...
			'orbit', 0, ...
			'ds', f.dx, 'dt', f.dx, ...
			'dfs', inf, ... % flat
			'dsd', 650, ...
			'dod', abs(igs.ny*igs.dy/2), ... % compressed on detector
			'down', f.down);

	A = Gtomosyn(cg, ig, 'chat', 0);
end

	x = ig.circ;
if 0, printm 'time'
	cpu etic
	y = A * x;
	cpu etoc 'A * x time:'
	% 160s for down=1 on jf27 (imac)
end

% todo: compare to analytical - oblique views look funny
