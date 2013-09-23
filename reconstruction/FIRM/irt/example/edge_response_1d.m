% edge_response_1d
% show how edge response is nonlinear for edge-preserving regularization

if ~isvar('ytrue'), printm 'setup'
	nx = 2^5;
	ny = nx;
	steps = [1 2 4 8];
	nstep = numel(steps);
	xtrue = zeros(nx, ny, 'single');
	ny_per = ny / nstep;
	for ii=1:nstep
		xtrue(nx/2+1:end, [1:ny_per] + (ii-1)*ny_per) = steps(ii);
	end

	psf = [1 2 1]';
	mask = true(nx,ny);
	A = Gblur(mask, 'psf', psf);
	ytrue = A * xtrue;

	if 1
		im('notick', xtrue, ' ')
		savefig edge_response_1d_xtrue
%	prompt
	end

	if 1
		im('notick', ytrue, ' ')
		savefig edge_response_1d_ytrue
%	prompt
	end
end

if ~isvar('R'), printm 'R'
	l2b = -0;
	bet = 2^l2b; % essentially product of beta * sigma^2
	f.pot = {'quad'};
	f.pot = {'hyper3', 1.0};

	R = Reg1(mask, 'beta', bet, 'type_penal', 'mat', ...
		'offsets', 1, 'pot_arg', f.pot); % 1d regularization
	if 1
		qpwls_psf(A, R, 1, mask);
	prompt
	end
end

if ~isvar('xh'), printm 'xh'
	xh = pwls_pcg1(xtrue(mask), A, 1, ytrue(:), R, 'niter', 30);
	xh = embed(xh, mask);

	if 1
		im('notick', xh, ' ')
	%	cbar, colormap(jet)
		savefig edge_response_1d_xh
	end
end

if 1 % plots
	tmp = zeros(nx, nstep);
	leg = cell(nstep,1);
	for ii=1:nstep
		tmp(:,ii) = xh(:, 1 + (ii-1)*ny_per) / steps(ii); % normalize
		leg{ii} = num2str(steps(ii));
	end
	plot(-nx/2:nx/2-1, tmp, '.-')
	legend(leg{:}, 2)
	axis([0*nx/2+[-1 1]*5 -0.1 1.1])
	ytick([0 1])
	xlabel 'horizontal location'
	ylabel 'normalized profile'
	savefig eps_c edge_response_1d_p1
end
