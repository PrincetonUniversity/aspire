  function fit = de_ftab_fit(sl, fm, varargin)
%|function fit = de_ftab_fit(sl, fm, [options])
%|
%| Fit to DE table fm(), suitable for subsequent interpolation / extrapolation.
%| Uses either classic polynomial basis functions,
%| or an experimental exponential model -log(sum_k p_k exp(-m_k . s))
%|
%| in
%|	sl	{L}		sample locations for each of L materials
%|	fm	[s1,...,sL,M]	DE tables f_m(s_1,...,s_L), size [(Ns),M]
%|
%| option
%|	'type'	'poly' | 'exp'	type of fit (default: 'poly')
%|	'wt'	{M}		fit weighting for each of M energies.
%|	'show'	1|0		plot?
%| options for poly
%|	'order'			polynomial order (default: 3)
%|	'maxdegree'		argument for poly_string() (default: [])
%| options for exp
%|	'kev'	[Ne,1]		kev *for fitting* (default: [10:5:200]')
%|	'mac'	[Ne,L]		*for fitting* (this or mtype are required)
%|	'mtype' {L}
%|	'macbar' [M,L]		if nonempty, use it to match derivative at s=0.
%|
%| out
%|	fit	strum		strum object for fitted f_m(s_1,...,s_L)
%|				fit.coef is [nbasis, M] for 'poly'
%|	methods:
%|	.fmfun(sll)		fm function evaluation, for stacked array sll
%|	.fgrad(sll)		fm gradient evaluation
%|	.fhess(sll)		fm hessian evaluation
%|	.curv2(sll, hfim)	curvatures of LS data-fit: (hfim - fm(s))^2/2
%|	.show_sp([en, sp])	plot true spectrum vs fitted spectrum
%|	.show_fm(sl, [fm])	mesh plot of fm and its fit
%|	.show_err(sl, fm)	mesh plot of fit error
%|	.mac_eff		effective mass atten coef based on fit
%|				(valid for 'exp' only)
%|
%| Copyright 2006-3-3, Jeff Fessler, University of Michigan

if nargin == 1 && streq(sl, 'test'), de_ftab_fit_test, return, end
if nargin < 2, help(mfilename), error(mfilename), end

arg.show = false;
arg.type = '';
arg.kev = [10:5:200]';
arg.mtype = {};
arg.mac = [];
arg.macbar = [];
arg.order = 3;
arg.maxdegree = [];
arg.wt = {};
arg.dc = false; % default: exclude dc (constant) from (polynomial) fit
arg = vararg_pair(arg, varargin);
if isempty(arg.type), arg.type = 'poly'; end

LL = length(sl);

if ndims(fm) == LL
	MM = 1;
elseif ndims(fm) == LL + 1
	MM = size(fm,ndims(fm)); % last dimension is number of spectra
else
	error 'invalid fm size'
end

for ll=1:LL % make sure arguments match table dimensions
	if length(sl{ll}) ~= size(fm, ll)
		fail('dim mismatch %d', ll)
	end
end

if LL > 2, warn 'todo: only L<=2 materials tested', end

switch arg.type
case 'poly'
	[fit fmfun fgrad fhess curv2] = de_ftab_fit_poly(sl, fm, MM, arg.wt, ...
		arg.order, arg.maxdegree, arg.dc);

case 'exp'
	fit = de_ftab_fit_exp(sl, fm, MM, arg.wt, ...
		arg.mtype, arg.mac, arg.kev, arg.macbar);
	fmfun = @de_ftab_fit_exp_eval;
	fgrad = @de_ftab_fit_exp_grad;
	fhess = @de_ftab_fit_exp_hess;
	curv2 = @de_ftab_fit_exp_curv2;
	fit.mac_eff = zeros(MM, LL); % [M,L] mac effective
	for mm=1:MM
		fit.mac_eff(mm,:) = fit.mac{mm}' * fit.coef{mm};
	end

otherwise
	fail('unknown fit type %s', arg.type)
end

fit.LL = LL;
fit.MM = MM;
fit.type = arg.type;
fit.cs = struct([]); % place holder for curvature structure

meth = {'fmfun', fmfun, '(sll)'; ...
	'fgrad', fgrad, '(sll)'; ...
	'fhess', fhess, '(sll)'; ...
	'curv2', curv2, '(sll, hfim)'; ...
	'show_err', @de_ftab_fit_show_err, '(sl, fm)'; ...
	'show_fm', @de_ftab_fit_show_fm, '(sl, [fm])'; ...
	'show_sp', @de_ftab_fit_show_sp, '([en, sp])'; ...
	};
fit = strum(fit, meth);

if arg.show && im
	fit.show_fm(sl, fm);
end

if LL == 1 && streq(arg.type, 'exp')
	fit.cs = de_ftab_fit_exp_curv2_setup(fit, 'show', arg.show);
end

end % de_ftab_fit()


%
% de_ftab_fit_poly()
%
function [fit, ffun, fgrad, fhess, curv2] = ...
	 de_ftab_fit_poly(sl, fm, MM, wt, order, maxdegree, dc)

sll = ndgrid_jf('mat', sl{:});

% for fitting, up weight the no-bone part
% since soft tissue is most prevalant
if isempty(wt)
	wt{1} = 1 ./ (stackpick(sll,2) + 1);
	wt{2} = wt{1};
end

fit.basis_order = order;
fit.basis_maxdegree = maxdegree;
fit.basis_dc = dc;

[fit.basis_str fit.basis_d1_str fit.basis_d2_str] = ...
	poly_string(order, 'maxdegree', maxdegree, 'dc', dc);
fit.basis_func = inline(fit.basis_str, 'x', 'y');
fit.basis_d1 = inline(['[ ' fit.basis_d1_str ' ]'], 'x', 'y');
fit.basis_d2 = inline(['[ ' fit.basis_d2_str ' ]'], 'x', 'y');
fit.nbasis = length(fit.basis_func(0,0));

ss1 = col(stackpick(sll,1));
ss2 = col(stackpick(sll,2));
tmp = fit.basis_func(ss1, ss2); % [*Ns,nbasis]
fit.coef = zeros(fit.nbasis, MM);
for mm=1:MM
	fit.coef(:,mm) = (diag_sp(wt{mm}(:)) * tmp) ...
		\ (wt{mm}(:) .* col(fm(:,:,mm)));
end

ffun = 'reshape([fit.basis_func(sl{1}(:), sl{2}(:)) * fit.coef(:,1); fit.basis_func(sl{1}(:), sl{2}(:)) * fit.coef(:,2)], [size(sl{1}) 2])';
ffun = inline(ffun, 'fit', 'sl');
ffun = @de_ftab_fit_poly_eval;
fgrad = @() error('not done');
fhess = @() error('not done');
curv2 = @() error('not done');

end % de_ftab_fit_poly()


%
% de_ftab_fit_poly_eval()
%
function fm = de_ftab_fit_poly_eval(fit, sll)

% if fit.MM ~= 2
s1 = stackpick(sll,1);
s2 = stackpick(sll,2);
fm = [	fit.basis_func(s1(:), s2(:)) * fit.coef(:,1);
	fit.basis_func(s1(:), s2(:)) * fit.coef(:,2)];
fm = reshape(fm, [size(s1) 2]);

end % de_ftab_fit_poly_eval()


%
% de_ftab_fit_exp()
% todo: more documentation
% todo: use a set of exponents, not "mac"
%
function fit = de_ftab_fit_exp(sl, fm, MM, wt, mtype, mac, kev, macbar)

sll = ndgrid_jf('mat', sl{:}); % [(Ns)]

if isempty(mac)
	if isempty(mtype), error 'mac or mtype required', end
	mac = xray_read_atten(mtype, kev); % [Ne,L]
else
	if size(mac,1) ~= length(kev), fail 'size mismatch', end
end

LL = length(sl);

if isempty(wt), wt = num2cell(ones(MM,1)); end

Ab = exp(-reshapee(sll, [], LL) * mac'); % [*Ns,Ne] "over-complete" basis

fit.kev = cell(1,MM);
fit.mac = cell(1,MM);
fit.coef = cell(1,MM);
for mm=1:MM
	if LL == 1 && ~isempty(macbar)
		warg = {'inprodv', mac' / macbar(mm,1)}; % constrain deriv. at 0
	else
		warg = {};
	end
	if MM == 1 % todo: kludgy
		dat = fm;
	else
		dat = stackpick(fm,mm); % [(Ns)]
	end
	y = exp(-dat);
	Wh = spdiag(sqrt(wt{mm}(:)), 'nowarn');
	% initial coefficients for each candidate energy
	x = wls_simplex(Ab, y(:), Wh, [], warg{:}); % [Ne,1]

	ie = x > 1e-6; % find key energies
	fit.kev{mm} = kev(ie);
	fit.mac{mm} = mac(ie,:); % [Ne,L] (now Ne may be smaller than before)

	if LL == 1 && ~isempty(macbar) % constrain derivative at 0
		warg = {'inprodv', fit.mac{mm}' / macbar(mm,1)};
	else
		warg = {};
	end

	A = exp(-reshapee(sll, [], LL) * fit.mac{mm}'); % [*Ns,Ne] final basis
	% final coefficients at key enerties:
	fit.coef{mm} = wls_simplex(A, y(:), Wh, [], warg{:}); % [Ne,1]
end

end % de_ftab_fit_exp()


%
% de_ftab_fit_exp_eval()
% evaluate 
% in
%	sll	[(Ns),L]	stackup of s1,s2,...,s_L
% out
%	f	[(Ns),M]	stackup of f1,f2,...,f_M
%
function f = de_ftab_fit_exp_eval(fit, sll)
Ns = size(sll); LL = Ns(end); Ns = Ns(1:end-1);
sll = reshapee(sll, [], LL); % [*Ns,L]
MM = fit.MM;

persistent warned
if ~isvar('warned') || isempty(warned)
	warned = 0;
end

f = zeros(prod(Ns),MM);
for mm=1:MM
	A = exp(-sll * fit.mac{mm}'); % [*Ns,Ne]
	if 0 % todo cut
		tmp = A;
		A = 1;
		mac = fit.mac{mm};
		for ll=1:LL
			sl = sll(:,ll);
			A = A .* exp(-sl * mac(:,ll)'); % [*Ns,Ne]
		end
		equivs(A, tmp)
	end
	tmp = -log(A * fit.coef{mm}); % [*Ns,1]
	f(:,mm) = tmp;

	if any(sll(:) < 0) % trick for negatives
		if LL == 1
			A = exp(-abs(sll) * fit.mac{mm}'); % [*Ns,Ne]
			tmp = -log(A * fit.coef{mm}); % [*Ns,1]
			f(:,mm) = -sign(sll) .* log(A * fit.coef{mm}); % [*Ns,1]
%			bad = sll < 0;
%			f(bad,mm) = sll(bad) * fit.mac_eff(mm,1);
		else
			if ~warned
				warned = 1;
				warn 'negative trick not done for L>1'
			end
		end
	end

end
f = reshape(f, [Ns MM]);

end % de_ftab_fit_exp_eval()


%
% de_ftab_fit_exp_grad()
% evaluate gradient of f for each of the given s vectors.
% in
%	sll	[(Ns),L]	stackup of s1,s2,...,s_L
% out
%	g	[(Ns),L,M]	stackup of gradients of f(s)
%
function g = de_ftab_fit_exp_grad(fit, sll)
Ns = size(sll); LL = Ns(end); Ns = Ns(1:end-1);
sll = reshape(sll, [], LL); % [*Ns,L]
MM = fit.MM;

persistent warned
if ~isvar('warned') || isempty(warned)
	warned = 0;
end

g = zeros(prod(Ns), LL, MM);
for mm=1:fit.MM
	alf = fit.coef{mm}; % [Ne,1]
	mac = fit.mac{mm}; % [Ne,L]
	Ne = length(alf); % # exponential terms in the fit

	if any(sll(:) < 0) % trick for negatives
		if LL == 1
			sll = abs(sll);
%			bad = sll < 0;
%			g(bad,1,mm) = fit.mac_eff(mm,1);
		else
			if ~warned
				warned = 1;
				warn 'negative trick not done for L>1'
			end
		end
	end

	A = exp(-sll * mac'); % [*Ns,Ne]
	Q = A .* repmat(alf', [prod(Ns) 1]); % [*Ns,Ne]
	denom = sum(Q,2); % [*Ns,1]
	Q = Q ./ repmat(denom, [1 Ne]);
	g(:,:,mm) = Q * mac;

end
g = reshape(g, [Ns LL MM]);

end % de_ftab_fit_exp_grad()


%
% de_ftab_fit_exp_hess()
% evaluate hessian of f for each of the given s vectors.
% in
%	sll	[(Ns),L]	stackup of s1,s2,...,s_L
% out
%	h	[(Ns),L,L,M]	stackup of hessians of f(s)
%
function h = de_ftab_fit_exp_hess(fit, sll)
Ns = size(sll); LL = Ns(end); Ns = Ns(1:end-1);
sll = reshape(sll, [], LL); % [*Ns,L]
MM = fit.MM;

persistent warned
if ~isvar('warned') || isempty(warned)
	warned = 0;
end

h = zeros(prod(Ns), LL, LL, MM);
for mm=1:fit.MM
	mac = fit.mac{mm}; % [Ne,L]
	alf = fit.coef{mm}; % [Ne,1]
	Ne = length(alf); % # exponential terms in the fit

	if any(sll(:) < 0) % trick for negatives
		if LL == 1
			sll_sign = sign(sll);
			sll = abs(sll);
%			bad = sll < 0;
%			hm(bad,1,1) = 0;
		else
			if ~warned
				warned = 1;
				warn 'negative trick not done for L>1'
			end
		end
	end

	A = exp(-sll * fit.mac{mm}'); % [*Ns,Ne]

	Q = A .* repmat(alf', [prod(Ns) 1]); % [*Ns,Ne]
	denom = sum(Q,2); % [*Ns,1]
	Q = Q ./ repmat(denom, [1 Ne]);
	g = Q * mac; % [*Ns,L] gradient

	% loop over LL because it is smaller than *Ns or Ne
	hm = zeros(prod(Ns), LL, LL);
	for l1=1:LL
		for l2=1:LL
			hm(:,l1,l2) = g(:,l1) .* g(:,l2) ...
				- Q * (mac(:,l1) .* mac(:,l2));
		end
	end

	if any(sll(:) < 0) % trick for negatives
		if LL == 1
			hm(:,1,1) = hm(:,1,1) .* sll_sign;
%			bad = sll < 0;
%			hm(bad,1,1) = 0;
		end
	end

	h(:,:,:,mm) = hm;

end
h = reshape(h, [Ns LL LL MM]);

end % de_ftab_fit_exp_hess()


%
% de_ftab_fit_exp_curv2_setup()
% tabulate maximum of (f'(s))^2 + (-d^/ds^2 f(s)) (hfim - f(s)),
% needed for surrogate curvature for LS cost function
%
function cs = de_ftab_fit_exp_curv2_setup(fit, varargin)

cs.show = 0;
cs.hfim_list = linspace(-5, 10, 31);
cs = vararg_pair(cs, varargin);
cs.show = cs.show & im;

LL = fit.LL;
if LL~=1, fail 'only LL=1 done', end
sl = linspace(-5, 5, 201)';
%sl = ndgrid_jf('mat', sl{:});
f = fit.fmfun(sl); % [*Ns,M]
g = fit.fgrad(sl);
g = squeeze(g); % [*Ns,M]
h = fit.fhess(sl);
h = squeeze(h); % [*Ns,M]

cs.ctable = zeros(length(cs.hfim_list), fit.MM);
for ii = 1:length(cs.hfim_list)
	hfim = cs.hfim_list(ii);
	tmp = g.^2 - h .* (hfim - f); % [*Ns,M]

	[cmax where] = max(tmp); % [1,M]
	if cs.show
		plot(sl, tmp, '-', sl(where), cmax, 'o')
	drawnow
%	prompt
	end
	if (any(where == 1) || any(where == length(sl))) && hfim > 0
		warn('cmax not in middle')
		keyboard
	end

	cs.ctable(ii,:) = cmax;
end
if cs.show, prompt, end

% curvature function: independent of sil because it is "max"
cs.curvf = @(cs, mm, sil, hfim) ...
	interp1(cs.hfim_list, cs.ctable(:,mm), hfim, 'linear', 'extrap');

if cs.show
	fh = linspace(-6, 12, 101);
	plot(cs.hfim_list, cs.ctable, 'o', fh, cs.curvf(cs, 1, 0, fh), '-')
	xlabel 'hat f im', ylabel 'curv2'
prompt
end

return

% old way:

tmp = abs(g.^2 + h .* f);
[cmax where] = max(tmp); % [1,M]
if cs.show
	plot(sl, tmp, '-', sl(where), cmax, 'o')
prompt
end
if any(where == 1) || any(where == length(sl))
	fail('cmax not in middle')
end

tmp = abs(h);
[vmax where] = max(tmp); % [1,M]
if any(where == 1) || any(where == length(sl))
	fail('vmax not in middle')
end
if cs.show
	plot(sl, tmp, '-', sl(where), vmax, 'o')
end

curv_max = {cmax, vmax};

end % de_ftab_fit_exp_curv2_setup()


%
% de_ftab_fit_exp_curv2()
% surrogate curvatures for LS costfunction for each of the given s vectors.
% in
%	sil	[(nd),L]	component material density integrals
%	hfim	[(nd),M]	\hat{f}_im (log estimates)
% out
%	curv	[(nd),M]	curvatures for LS cost function
%
function curv = de_ftab_fit_exp_curv2(fit, sil, hfim)
LL = fit.LL;
MM = fit.MM;
nd = size(hfim); if MM ~= nd(end), error 'hfim size', end; nd = nd(1:end-1);
sil = reshapee(sil, [], LL); % [*nd,M]
hfim = reshapee(hfim, [], MM); % [*nd,M]
if size(sil,1) ~= prod(nd), error 'sil size', end
curv = zeros(prod(nd), MM);

cs = fit.cs;

%curv_max = fit.curv_max;
%cmax = curv_max{1}; % [1,M]
%vmax = curv_max{2}; % [1,M]
for mm=1:fit.MM
%	mac = fit.mac{mm}; % [Ne,L]
%	vmax = (max(mac) - min(mac))^2 / 4;
%	curv(:,mm) = max(mac)^2 + vmax * hfim(:,mm) + cmax(mm);
%	curv(:,mm) = vmax(mm) * abs(hfim(:,mm)) + cmax(mm);
	curv(:,mm) = cs.curvf(cs, mm, sil, hfim(:,mm));
end
curv = reshape(curv, [nd MM]);

end % de_ftab_fit_exp_curv2()


%
% de_ftab_fit_show_sp()
% compare true spectra to fitted spectra
%
function out = de_ftab_fit_show_sp(fit, en, sp)
if nargin < 1, fail 'de_ftab_fit_show_sp(fit, en, sp)', end

if ~streq(fit.type, 'exp'), printm 'show_sp only done for exp', return, end
if im
	clf, pl = (fit.MM+1)*100 + 10 + 1;
	if nargin > 1
		subplot(pl)
		plot(en, sp * diag(1 ./ max(sp)))
		tmp = num2cell(char([1:fit.MM] + '0'));
		legend(tmp{:})
		axisx(minmax(en))
	end

	if isfield(fit, 'kev')
		for mm=1:fit.MM
			subplot(pl+mm)
			bar(fit.kev{mm}, fit.coef{mm})
			xtick(fit.kev{mm}(fit.coef{mm} > 0))
			axis tight
			if nargin > 1, axisx(minmax(en)), end
		end
	else
		warn 'kev unknown'
	end
end

if nargout, out = []; end

end % de_ftab_fit_show_sp()


%
% de_ftab_fit_show_fm()
% compare fit to sampled fm
%
function out = de_ftab_fit_show_fm(fit, sl, fm, sll)
if nargin < 2, error 'de_ftab_fit_show_fm(fit, sl, fm)', end

do_fm = nargin > 2;
if ~isvar('sll') || isempty(sll)
	sll = ndgrid_jf('mat', sl{:});
end

fh = fit.fmfun(sll);

switch length(sl) % LL
case 1
	s1 = sl{1};
	im clf
	arg = {};
	leg = {};
	order = 'bgrcm';
	if do_fm
		for mm=1:fit.MM
			arg = {arg{:}, s1, fm(:,mm), [order(mm) '.']};
			leg = {leg{:}, sprintf('true m=%d', mm)};
		end
	end
	for mm=1:fit.MM
		arg = {arg{:}, sll, fh(:,mm), [order(mm) '-']};
		leg = {leg{:}, sprintf('fit m=%d', mm)};
	end
	plot(arg{:})
	xlabel 's1'
	ylabel 'fm(s1)'
	legend(leg{:}, 4)

case 2
	s1 = sl{1};
	s2 = sl{2};
	smax(1) = max(s1);
	smax(2) = max(s2);
	if do_fm
		fmax = max(fm(:));
	else
		fmax = max(fh(:));
	end

	ax = [0 smax(1) 0 smax(2) 0 fmax];
	cax = [0 fmax];

	im clf
	im('pl', 2, fit.MM)
	for mm=1:fit.MM
		if do_fm
		show(mm, s1, s2, fm(:,:,mm), ax, cax, sprintf('f_%d(s)', mm))
		end
		m0 = mm + fit.MM;
		show(m0, s1, s2, fh(:,:,mm), ax, cax, sprintf('f_%d fit', mm))
	end
%	% text(-60, 10, '[cm^2/g]')
otherwise
	fail 'not done'
end

if nargout, out = []; end

end % de_ftab_fit_show_fm()


%
% de_ftab_fit_show_err()
% show fit errors, and relate to HU
% f = mac s, mac(70kev,H2O) = 0.2 cm^2 / g = 1000 HU
% mh = mac * 1000 HU / (0.2 g/cm^2)
% so Df/Dmh = Df/Dmac * Dmac/Dmh = (50 g/cm^2) (0.2 cm^2/g) / 1000 HU = 1/100 HU
% so Dmh = 100 HU * Df for 50 cm of H2O.
%
function out = de_ftab_fit_show_err(fit, sl, fm)
if nargin < 3, error 'de_ftab_fit_show_err(fit, sl, fm)', end

sll = ndgrid_jf('mat', sl{:});

fh = fit.fmfun(sll);
err = fh - fm;
printm('worst model error = %g of %g', max(abs(err(:))), max(fm(:)))
printm('worst model error = %g HU over 50cm H2O', 100*max(abs(err(:))))
for mm=1:fit.MM
	ee = stackpick(err,mm);
	printm('worst error (mm=%d): %g', mm, max(abs(ee(:))))
end

if max(abs(err(:))) == 0, return, end
if 0
	e1 = err(:,:,1);
	e2 = err(:,:,2);
	disp([minmax(e1); minmax(e2)]')
	printm('worst error1 %g', max(col(abs(err(:,1,:)))))
	printm('worst error2 %g', max(col(abs(err(:,2,:)))))
end
%err = abs(err);

switch length(sl) % LL
case 1
	s1 = sl{1};
	arg = {};
	leg = {};
	order = 'bgrcm';
	for mm=1:fit.MM
		arg = {arg{:}, s1, err(:,mm), [order(mm) '-']};
		leg = {leg{:}, sprintf('m=%d', mm)};
	end
	plot(arg{:})
	xlabel 's1'
	ylabel 'error'
	legend(leg{:}, 4)

case 2
	s1 = sl{1};
	s2 = sl{2};

	elim = minmax(err(:))';
	%elim = [-1 1] * 0.01; % +/- 1 HU
	ax = [0 max(s1) 0 max(s2) elim];
	im clf, im('pl', 1, fit.MM)
	for mm=1:fit.MM
		show(mm, s1, s2, err(:,:,mm), ax, elim, sprintf('f_%d error', mm))
	end
otherwise
	fail 'not done'
end

if nargout, out = []; end

end % de_ftab_fit_show_err()


%
% show()
%
function show(pl, x, y, f, ax, cax, ti)
if ~im, return, end
im('subplot', pl)
if 1
	mesh(x,y,f')
	colormap hsv, caxis(cax), cbar
	axis(ax)
	xtick, ytick, ztick, zwhite
else
	im(x,y,f), cbar
	xtick, ytick
end
xlabel 's_1', ylabel 's_2', title(ti)

end % show()


%
% de_ftab_fit_test()
%
function de_ftab_fit_test

%stype = 'mono,70';
stype = 'ps1';
stype = 'poly1,60';
stype = 'poly1,60,100,160'; % stress test
xrs = xray_read_spectra(stype);

list.sl{1} = linspace(0, 50, 26); % coarse for fast fitting
list.sl{2} = linspace(0, 30, 31);
list.slfine{1} = linspace(0, 5, 261); % fine for evaluating
list.slfine{2} = linspace(0, 3, 511);
list.mtype = {'water', 'bone'};

for LL=1:2
	sl = {list.sl{1:LL}};
	mtype = {list.mtype{1:LL}};

	mas = xray_read_mac(mtype);
	mac = mas.mac(xrs.en);
	if im
		clf, semilogy(xrs.en, mac), legend(mtype{:})
	end
	sll = ndgrid_jf('mat', sl{:});
	fm = de_ftab_fm(sll, mac, xrs.Ide);
	fit = de_ftab_fit(sl, fm, 'type', 'exp', 'mtype', mtype, 'show', 1);
%	fit = de_ftab_fit(sl, fm, 'type', 'poly')

	fh = fit.fmfun(sll);
	g = fit.fgrad(sll);
	h = fit.fhess(sll);

	if 0 && LL == 2 % test gradient
		slfine = {list.slfine{1:LL}};
		slf = ndgrid_jf('mat', slfine{:});
		fh = fit.fmfun(slf);
		g = fit.fgrad(slf);
		i1 = 5;
		i2 = 9;
		d1 = slfine{1}(i1+1) - slfine{1}(i1);
		d2 = slfine{2}(i2+1) - slfine{2}(i2);
		gh(:,1) = (fh(i1+1,i2,:) - fh(i1,i2,:)) / d1;
		gh(:,2) = (fh(i1,i2+1,:) - fh(i1,i2,:)) / d2;
		pr transpose(gh)
		gg = squeeze(g(i1,i2,:,:));
		pr gg
		pr transpose(fit.mac_eff) % at origin
%		gg - fit.mac_eff'
	end

	if 1 && LL == 1 % test surrogate and its curvature
		s1 = 5;
		f1 = fit.fmfun(s1);
		g1 = fit.fgrad(s1);
		g1 = squeeze(g1)'
		hf1 = f1 * (1.4);
		curv = fit.curv2(s1, hf1)
		if 0 % test newton curvature - hugs closely but crosses
			h1 = fit.fhess(s1);
			h1 = squeeze(h1)';
			curv = g1.^2 - h1 .* (hf1 - f1)
		end
		s = linspace(-8,8,201)';

		if im
			fit.show_fm(sl, fm, s); axisx(minmax(s))
		end

		f = fit.fmfun(s);
		for mm=1:fit.MM
			q1 = (hf1(mm) - f1(mm))^2/2 ...
				+ g1(mm) * (f1(mm)-hf1(mm)) * (s-s1) ...
				+ curv(mm)/2 * (s-s1).^2;
			cost = (hf1(mm)-f(:,mm)).^2/2;
			if im
				plot(s, cost, '-', s, q1, '--')
				axis([minmax(s)' minmax(cost)'])
			end
		prompt
		end

%		keyboard
	prompt
	end

	if im
		fit.show_sp(xrs.en, xrs.sp);
		prompt
		fit.show_fm(sl, fm);
		prompt
		fit.show_err(sl, fm);
		prompt
	end

end

end % de_ftab_fit_test()
