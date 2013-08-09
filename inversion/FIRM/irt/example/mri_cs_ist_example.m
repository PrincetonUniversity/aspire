% mri_cs_ist_example.m
% Compressed-sensing for MRI via iterative soft thresholding (IST)
% interactive demo with slider to select soft threshold value
% Copyright 2010-02-23, Jeff Fessler, University of Michigan

if ~isvar('A'), printm 'setup Gdft object'
	Nd = [1 3/2] * 2^7;
	rand('state', 1)
	samp = rand(Nd) > 0.3;
	mask = true(Nd);
%	mask(1) = false; % todo: stress
	A = Gdft('mask', mask, 'samp', samp);

	im clf, im(samp)
prompt
end

if ~isvar('U'), printm 'setup Gwave object'
	U = Gwave1('mask', mask, 'noCA', 0, 'level', 3);
%	U = Gwave1('mask', true([2^4 2^3]), 'noCA', 0, 'level', 3); im(U')
	U = U';
%	pr size(U)
	ic = 3500;
%	ic = randint(1,1,[1 size(U,2)])
	ic = min(ic,prod(Nd));
	im(embed(U(:,ic), mask))
prompt
end

clim = [0 9];
if 0 | ~isvar('xt'), printm 'setup object'
	ig = image_geom('nx', Nd(1), 'ny', Nd(2), 'dx', 1);
	xt = ellipse_im(ig, 'shepplogan-emis', 'oversample', 2);
	im(abs(xt), '|x| true', clim), cbar
prompt

	% show wavelet coefficients (cheat: helps determine scale)
	im(reshape(U' * xt, Nd), [-15 15]), cbar
prompt
end

	soft = @(t,a) (t - a * sign(t)) .* (abs(t) > a);
%	soft = @(z,a) z .* (t - a * sign(t)) .* (abs(t) > a); % todo!
	thresh = 2^-4;

if 0 % explore
	tmp = U * soft(U' * xt, thresh);
	tmp = embed(tmp, mask);
	im([xt; tmp; tmp-xt], clim)
return
end

if ~isvar('yi'), printm 'setup data'
	yb = A * xt(mask);
	sig = 30;
	yi = yb + sig * randn(size(yb));
	snr = 20 * log( norm(yb(:)) / norm(yi(:) - yb(:)) );
	pr snr
	im(log(max(abs(fftshift(embed(yi,samp))),1)), [1 10]), cbar
prompt
end


if ~isvar('xf'), printm 'basic inverse FFT reconstruction'
	xf = A' * yi / prod(Nd);
	xf = embed(xf, mask);
	im(abs(xf), '|x| "zero-pad"', clim), cbar
prompt
end

if ~im, return, end

if 1 || ~isvar('xi'), printm 'IST reconstruction'

	im plc 2 1, im(1, abs(xf), clim)
	hu = uicontrol('style', 'slider', 'string', 'value', ...
		'units', 'normalized', 'position', [0.1 0.0 0.8 0.03], ...
		'value', 0);

	hb = uicontrol('style', 'togglebutton', 'string', 'stop',...
		'units', 'normalized', 'position', [0.01 0.0 0.08 0.03]);
	drawnow
prompt

	tfun = @(x) 2^-8 + 2^-1 * x;
	thresh_save = nan;

	% todo: somehow threshold only the detail coefficients
	curv = prod(Nd); % because standard DFT

%	reg = thresh * curv;
	xi = xf(mask);
	while(1)
		if (get(hb, 'value')) % stop button
			break
		end
		thresh = tfun( get( hu, 'value') );
		if thresh_save ~= thresh
			thresh_save = thresh;
			iter = 1;
		end
		tmp = xi + 1/curv * A' * (yi - A * xi);
		xi = U * soft(U' * tmp, thresh);
		if ~rem(iter, 4)
			% thresholded truth for comparison
			tmp = embed(U * soft(U' * xt, thresh), mask);
			tmp = [tmp; embed(abs(xi), mask)];
			im(1, tmp, clim, ' '), cbar
			fun = @(p,s,v) text(p*ig.nx/2, -4, sprintf(s, v), 'horiz', 'cent');
			fun(1, 'True: threshold = %g', thresh)
			fun(3, 'IST %4d', iter)

			wl = [reshape(U' * xt(:), Nd); reshape(U' * xi(:), Nd)];
			im(2, abs(wl), [0 4])
			drawnow
		end
		iter = iter + 1;
	end
	xi = embed(xi, mask);
end

if 0
	im plc 1 3
	im(1, xt)
	im(2, xf)
	im(3, xi)
end
