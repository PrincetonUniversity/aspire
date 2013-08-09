 function psino = rebin_fan2par(fsino, sf, sp, varargin)
%function psino = rebin_fan2par(fsino, sf, sp, varargin)
%
% Rebin fan-beam sinogram into parallel-beam (or mojette) sinogram.
% (Also useful for parallel-to-mojette rebinning.)
%
% in
%	fsino	[ns,nbeta]	fan-beam sinogram (or possibly parallel)
%	sf			sino_geom() for fan-beam input
%	sp			sino_geom() for parallel-beam output
%
% options
%	'ob'			set to 1 to create (Fatrix) object
%	's_interp'	default: {'order', 3, 'ending', 'zero'}
%	'beta_interp'	default: {'order', 3, 'ending', 'periodic'}
%
% out
%	psino	[nr,nphi]	parallel-beam sinogram (or possibly mojette)
%
% Copyright 2005-12-10, Jeff Fessler, The University of Michigan

if nargin == 1 && streq(fsino, 'test'), rebin_fan2par_test, return, end
if nargin < 3, help(mfilename), error(mfilename), end

% defaults
arg.ob = false;
arg.s_interp = {'order', 3, 'ending', 'zero'};
arg.beta_interp = {'order', 3, 'ending', 'periodic'};
arg = vararg_pair(arg, varargin);

if isempty(sf.nb), sf.nb = size(fsino,1); end
if isempty(sf.na), sf.na = size(fsino,2); end
arg.dimi = sf.dim;
arg.dimo = sp.dim;

if streq(sp.type, 'par')
	is_mojette = 0;
	dr = sp.dr;
elseif streq(sp.type, 'moj')
	is_mojette = 1;
	dr = sp.dx; % trick:
else
	error 'output sino must be "par" or "moj"'
end

arg.sf = sf;

if streq(sf.type, 'fan') % fan->(par|moj)
	[arg.r_ob arg.phi_ob arg.flag180] = rebin_fan2par_setup(...
		sf.nb, sf.ds, sf.offset_s, ...
		sf.na, sf.orbit_start, sf.orbit, ...
		sf.dso, sf.dsd, sf.dfs, ...
		sp.nb, dr, sp.offset_r, ...
		sp.na, sp.orbit_start, sp.orbit, ...
		is_mojette, arg.s_interp, arg.beta_interp);

elseif streq(sf.type, 'par') % par->(moj|par)
	[arg.r_ob arg.phi_ob arg.flag180] = rebin_fan2par_setup(...
		sf.nb, sf.dr, sf.offset_r, ...
		sf.na, sf.orbit_start, sf.orbit, ...
		inf, inf, 0, ...
		sp.nb, dr, sp.offset_r, ...
		sp.na, sp.orbit_start, sp.orbit, ...
		is_mojette, arg.s_interp, arg.beta_interp);

else
	error([sp.type ' to ' sf.type ' not done'])
end

if arg.ob
	psino = rebin_fan2par_ob(arg);
else
	psino = rebin_fan2par_arg(arg, fsino);
end


%
% rebin_fan2par_ob()
%
function ob = rebin_fan2par_ob(arg)

dim = [prod(arg.dimo) prod(arg.dimi)];
ob = Fatrix(dim, arg, 'caller', 'rebin_fan2par', ...
        'forw', @rebin_fan2par_arg, 'back', @error);


%
% rebin_fan2par_arg()
%
function sino = rebin_fan2par_arg(arg, sino)

if size(sino,1) == prod(arg.dimi)
	sino = reshape(sino, arg.dimi(1), arg.dimi(2), []);
	flag_column = 1;
else
	flag_column = 0;
end

sino = (arg.phi_ob * sino.').';

if arg.flag180
	sino = rebin_fan2par_inlace(arg, sino);
end

sino = arg.r_ob * sino;
if flag_column
	sino = reshape(sino, prod(arg.dimo), []);
end


%
% rebin_fan2par_inlace()
% interlace opposing views (for quarter-detector offset)
%
function sino = rebin_fan2par_inlace(arg, sino)

if ndims(sino) > 2, error 'multisino not done', end
ns = arg.dimi(1);
nphi = arg.dimo(2);
t1 = sino(:,1:nphi); % [ns,nphi]
t2 = flipdim(sino(:,nphi + [1:nphi]), 1); % [ns,nphi]
if arg.sf.offset == 1.25 % trick
	t1 = [t1 ; t2([ end-1 end],:)]; % [ns+2,nphi]
	t2 = [t1([1 2],:) ; t2]; % [ns+2,nphi]
	ns = ns + 2;
elseif arg.sf.offset ~= 0.25
	error 'bug'
end
sino = reshape([t1(:)'; t2(:)'], 2*ns, nphi, []);


%
% rebin_fan2par_setup()
%
function [r_ob phi_ob flag180] = rebin_fan2par_setup(...
	ns, ds, offset_s, ...
	nbeta, beta_start, beta_orbit, ...
	dso, dsd, dfs, ...
	nr, dr, offset_r, ...
	nphi, phi_start, phi_orbit, ...
	is_mojette, ... % 1 output sinogram is to be mojette
	s_interp, beta_interp)

if dfs, error 'flat fan not done', end

if phi_orbit == 180 && beta_orbit == 360
	flag180 = 1; % trick: handle 180 -> 360
else
	flag180 = 0;
end

phi_start = deg2rad(phi_start);
phi_orbit = deg2rad(phi_orbit);
beta_start = deg2rad(beta_start);
beta_orbit = deg2rad(beta_orbit);
phi = phi_start + phi_orbit * [0:nphi-1]' / nphi;
%bet = beta_start + beta_orbit * [0:nbeta-1]' / nbeta;

% angular interpolator
ws = (ns-1)/2 + offset_s;

if isinf(dsd) % parallel - no need for angular interp if orbits match
	if phi_start == beta_start && phi_orbit == beta_orbit && nphi == nbeta
		phi_ob = 1;
	else
		error 'parallel beam angular interpolation not done'
	end

else % fan

	if flag180 % expand desired phi's from [x,x+180] to [x,x+360]
		phi2 = phi_start + 2 * phi_orbit * [0:(2*nphi-1)]' / (2*nphi);
	else
		if phi_orbit ~= deg2rad(360) || beta_orbit ~= deg2rad(360)
			error 'todo: only 360 done - ask jeff'
		end
		phi2 = phi;
	end

	s = ([0:ns-1]' - ws) * ds;

	bet = outer_sum(phi2, -s / dsd); % beta = phi - gam
	bet_int = nbeta / beta_orbit * (bet - beta_start);

	phi_ob = bspline_1d_interp(zeros(nbeta,ns), bet_int, ...
		beta_interp{:}, 'ob', 1);
end

% radial interpolator

wr = (nr-1)/2 + offset_r;
if is_mojette
	dr = dr * max(abs(cos(phi)), abs(sin(phi)))';
end
r = ([0:nr-1]' - wr) * dr; % trick: [nr,1] or [nr,nphi]
if isinf(dsd)
	s = r;
else
	s = dsd * asin(r / dso);
end

if flag180 % trick: due to interlacing, the effective "s" sampling changes
	if offset_s ~= 0.25 && offset_s ~= 1.25
		error 'only 0.25 and 1.25 implemented now; but generalizable'
	end
	ns = ns * 2 + 4 * (offset_s - 0.25);
	offset_s = 0;
	ws = (ns-1)/2 + offset_s;
	ds = ds / 2;
end
s_int = s / ds + ws;
r_ob = bspline_1d_interp(zeros(nr, nphi), s_int, s_interp{:}, 'ob', 1);


%
% test both rebin_fan2par() and par2fan_rebin()
%
function rebin_fan2par_test

down = 4;
dx = down/2;

gp = sino_geom('par', 'nb', 1096/down, 'na', 800/down, ...
	'dr', 0.5*down, 'offset_r', 0, 'orbit', 180);
gm = sino_geom('moj', 'nb', gp.nb, 'na', gp.na, ...
	'dx', dx, 'offset_r', 0, 'orbit', 180);
gf = sino_geom('fan', 'ns', 888/down, 'nbeta', 984/down, ...
	'dsd', 949, 'dod', 408, ...
	'ds', down, 'offset_s', 1.25); % quarter detector

% analytical sinograms
ell = [0 20 180 150 0 1];
oversample = 4;
par = ellipse_sino(gp, ell, 'oversample', oversample);
moj = ellipse_sino(gm, ell, 'oversample', oversample);
fan = ellipse_sino(gf, ell, 'oversample', oversample);

% sanity check for par->par
if 0
	g0 = sino_geom('fan', 'nb', gp.nb, 'na', gp.na, ...
		'dsd', inf, 'dod', 0, 'orbit', gp.orbit, ...
		'ds', gp.dr, 'offset_s', gp.offset_r);
	par2par = rebin_par2fan(par, gp, g0);
	max_percent_diff(par, par2par)
	par2par = rebin_fan2par(par, g0, gp);
	max_percent_diff(par, par2par)
return
end

if 1 % test rebin (normal)
	cpu etic
	fan2par = rebin_fan2par(fan, gf, gp);
	cpu etoc 'fan2par time'

	cpu etic
	par2fan = rebin_par2fan(par, gp, gf);
	cpu etoc 'par2fan time'

	max_percent_diff(par, fan2par)
	max_percent_diff(fan, par2fan)
	im clf, im pl 2 3
	im(1, par, 'par'), cbar
	im(4, fan, 'fan'), cbar
	im(2, fan2par, 'fan2par'), cbar
	im(5, par2fan, 'par2fan'), cbar
	im(3, fan2par-par, 'error'), cbar
	im(6, par2fan-fan, 'error'), cbar
prompt
end

% test rebin (mojette)
cpu etic
fan2moj = rebin_fan2par(fan, gf, gm);
cpu etoc 'fan2moj time'

cpu etic
moj2fan = rebin_par2fan(moj, gm, gf);
cpu etoc 'moj2fan time'

max_percent_diff(fan, moj2fan)
max_percent_diff(moj, fan2moj)

im(1, moj, 'moj'), cbar
im(4, fan, 'fan'), cbar
im(2, fan2moj, 'fan2moj'), cbar
im(5, moj2fan, 'moj2fan'), cbar
im(3, fan2moj-moj, 'error'), cbar
im(6, moj2fan-fan, 'error'), cbar
