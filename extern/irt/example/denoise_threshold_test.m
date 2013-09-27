%| denoise_threshold_test.m
%| Test denoising based on l_1 type penalty, cf thresholding
%| min_x 1/2 |y - x|^2 + \beta pot(x)
%|
%| Copyright 2005-4-22, Jeff Fessler, University of Michigan

yi = linspace(-10,10,101)';
mask = true(size(yi));
A = Gdiag(ones(size(mask))); % identity "matrix"

if 1
	f.l2b_q = 1;
%	f.l2b_n = 5;
	f.l2b_n = 1;
	f.cut = 5; % cutoff point
%	f.delta = 0.2;
	f.delta = f.cut / (1 + 2^f.l2b_n); % for broken parabola
	f.tik = (1 + 2^f.l2b_n) * f.delta;
	f.niter = 40;
%	f.type = 'cauchy';
	f.type = 'broken';

 if 0
	Rq = Robject(mask, 'type_denom', 'matlab', ...
		'offsets', 0, 'potential', 'quad', 'beta', 2^f.l2b_q);
 else
	Rq = Reg1(mask, 'type_denom', 'matlab', ...
		'offsets', 0, ... % trick for identity
		'beta', 2^f.l2b_q);
%		'pot_arg', {'quad'}, 'beta', 2^f.l2b_q);
 end

	xq = pwls_sps_os(0*yi(:), yi(:), [], A, Rq, ...
			2, [-inf inf], [], [], 1);
end

if 1
 if 0
	Rc = Robject(mask, 'type_denom', 'matlab', ...
		'offsets', 0, ... % trick for identity
		'potential', f.type, 'delta', f.delta, 'beta', 2^f.l2b_n);
 else
	Rc = Reg1(mask, 'type_denom', 'matlab', ...
		'offsets', 0, ... % trick for identity
		'pot_arg', {f.type, f.delta}, 'beta', 2^f.l2b_n);
 end

	xc = pwls_sps_os(0*yi(:), yi(:), [], A, Rc, ...
		f.niter, [-inf inf], [], [], 1);

 if 0
	Rh = Robject(mask, 'type_denom', 'matlab', ...
		'offsets', 0, ... % trick for identity
		'potential', 'hyper3', 'delta', f.delta, 'beta', 2^f.l2b_n);
 else
	Rh = Reg1(mask, 'type_denom', 'matlab', ...
		'offsets', 0, ... % trick for identity
		'pot_arg', {'hyper3', f.delta}, 'beta', 2^f.l2b_n);
 end

	xh = pwls_sps_os(yi(:), yi(:), [], A, Rh, ...
		f.niter, [-inf inf], [], [], 1);
end

if im
	clf, plot(yi, yi, ':', ...
		yi, xh(:,end), '-', ...
		yi, xc(:,end), '-.', ...
		yi, xq(:,end), '--')
%	axis equal
	axis square
	tik = [min(yi), -f.tik 0 f.tik max(yi)];
%%%%%%%%	xtick(tik), ytick(tik)
	grid
	legend('I', 'hyper', f.type, 'quad', 2)
end
