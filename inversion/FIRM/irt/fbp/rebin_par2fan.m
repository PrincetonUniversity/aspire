 function fsino = rebin_par2fan(psino, pgeom, fgeom, varargin)
%function fsino = rebin_par2fan(psino, pgeom, fgeom, varargin)
%
% Rebin parallel-beam (or mojette) sinogram into fan-beam sinogram.
% Also useful for mojette-to-parallel rebinning.
%
% in
%	psino	[nr,nphi]	parallel-beam sinogram
%	pgeom			sino_geom() for parallel-beam input
%	fgeom			sino_geom() for fan-beam output
%
% options
%	'ob'			set to 1 to create (Fatrix) object
%	'r_interp'	default: {'order', 3, 'ending', 'zero'}
%	'phi_interp'	default: {'order', 3, 'ending', 'periodic'}
%
% out
%	fsino	[ns,nbeta]	fan-beam sinogram (or possibly mojette)
%
% Copyright 2005-12-7, Jeff Fessler, The University of Michigan

if nargin == 1 && streq(psino, 'test'), rebin_par2fan_test, return, end
if nargin < 3, help(mfilename), error(mfilename), end

% defaults
arg.ob = false; % set to 1 to create (Fatrix) object
arg.r_interp = {'order', 3, 'ending', 'zero'};
arg.phi_interp = {'order', 3, 'ending', 'periodic'};

arg = vararg_pair(arg, varargin);

if isempty(pgeom.nb), pgeom.nb = size(psino,1); end
if isempty(pgeom.na), pgeom.na = size(psino,2); end
arg.dimi = [pgeom.nb pgeom.na];
arg.dimo = [fgeom.nb fgeom.na];

if streq(pgeom.type, 'par')
	is_mojette = 0;
	dr = pgeom.dr;
elseif streq(pgeom.type, 'moj')
	is_mojette = 1;
	dr = pgeom.dx;
else
	error 'input sino must be "par" or "moj"'
end

if streq(fgeom.type, 'fan') % (par|moj)->fan
	[arg.r_ob arg.phi_ob arg.flag180] = rebin_par2fan_setup(is_mojette, ...
		pgeom.nb, dr, pgeom.offset_r, ...
		pgeom.na, pgeom.orbit_start, pgeom.orbit, ...
		fgeom.nb, fgeom.ds, fgeom.offset_s, ...
		fgeom.na, fgeom.orbit_start, fgeom.orbit, ...
		fgeom.dso, fgeom.dsd, fgeom.dfs, ...
		arg.r_interp, arg.phi_interp);

elseif is_mojette && streq(fgeom.type, 'par') % moj->par
	[arg.r_ob arg.phi_ob arg.flag180] = rebin_par2fan_setup(true, ...
		pgeom.nb, dr, pgeom.offset_r, ...
		pgeom.na, pgeom.orbit_start, pgeom.orbit, ...
		fgeom.nb, fgeom.dr, fgeom.offset_r, ...
		fgeom.na, fgeom.orbit_start, fgeom.orbit, ...
		inf, inf, 0, ...
		arg.r_interp, arg.phi_interp);

else
	error([pgeom.type ' to ' fgeom.type ' not done'])
end

if arg.ob
	fsino = rebin_par2fan_ob(arg);
else
	fsino = rebin_par2fan_arg(arg, psino);
end


%
% rebin_par2fan_ob()
%
function ob = rebin_par2fan_ob(arg)

dim = [prod(arg.dimo) prod(arg.dimi)];
ob = Fatrix(dim, arg, 'caller', 'rebin_par2fan', ...
	'forw', @rebin_par2fan_arg, 'back', @rebin_par2fan_adj);


%
% rebin_par2fan_arg()
%
function sino = rebin_par2fan_arg(arg, sino)

if size(sino,1) == prod(arg.dimi)
	sino = reshape(sino, arg.dimi(1), arg.dimi(2), []);
	flag_column = 1;
else
	flag_column = 0;
end

if arg.flag180
	sino = [sino, flipdim(sino,1)];
end
sino = arg.r_ob * sino;
sino = (arg.phi_ob * sino.').';
if flag_column
	sino = reshape(sino, prod(arg.dimo), []);
end


%
% rebin_par2fan_setup()
% set up the objects needed for (repeated) interpolation
%
function [r_ob phi_ob flag180] = rebin_par2fan_setup(...
	is_mojette, ... % input sinogram mojette? trick: if so, dr is really dx
	nr, dr, offset_r, ...
	nphi, phi_start, phi_orbit, ...
	ns, ds, offset_s, ...
	nbeta, beta_start, beta_orbit, ...
	dso, dsd, dfs, ...
	r_interp, phi_interp)

if dfs, error 'flat fan not done', end

if phi_orbit == 180 && beta_orbit == 360 % trick: handle 180 -> 360
	if offset_r ~= 0
		error 'par/moj 180 -> fan 360 works only if offset_r=0'
	end
	flag180 = 1;
	phi_orbit = 2 * phi_orbit;
	nphi = 2 * nphi;
else
	flag180 = 0;
end
phi_start = deg2rad(phi_start);
phi_orbit = deg2rad(phi_orbit);
beta_start = deg2rad(beta_start);
beta_orbit = deg2rad(beta_orbit);
phi = phi_start + phi_orbit * [0:nphi-1] / nphi;
bet = beta_start + beta_orbit * [0:nbeta-1]' / nbeta;

% radial interpolation args

wr = (nr-1)/2 + offset_r;
ws = (ns-1)/2 + offset_s;
s = ([0:ns-1]' - ws) * ds;
if isinf(dsd)
	r = s;
	if ~is_mojette
		warning 'par2fan with parallel but not mojette?'
	end
else
	r = dso * sin(s / dsd);
end
if is_mojette
	dr = dr * max(abs(cos(phi)), abs(sin(phi)));
end
r_int = r * (1 ./ dr) + wr; % trick: [nr,1] or [nr,nphi]
r_ob = bspline_1d_interp(zeros(nr,nphi), r_int, r_interp{:}, 'ob', 1);

% angular interpolation args

if isinf(dsd) ...
	&& phi_start == beta_start && phi_orbit == beta_orbit && nphi == nbeta
	phi_ob = 1; % no angular interp if orbits match in parallel-beam!
return
end

if phi_orbit ~= 2*pi || beta_orbit ~= 2*pi
	error 'only 360 done - ask jeff'
end

if isinf(dsd)
	phi = bet;
%	if beta_start >= phi_start ...
%		&& beta_start + beta_orbit <= phi_start + phi_orbit ...
%		&& beta_orbit >= 0 && phi_orbit >= 0
else
	phi = outer_sum(bet, s / dsd);
end
phi_int = nphi / phi_orbit * (phi - phi_start);

phi_ob = bspline_1d_interp(zeros(nphi,nr), phi_int, phi_interp{:}, 'ob', 1);


%
% rebin_par2fan_adj()
%
function sino = rebin_par2fan_adj(arg, sino)

if size(sino,1) == prod(arg.dimo)
	sino = reshape(sino, arg.dimo(1), arg.dimo(2), []);
	flag_column = 1;
else
	flag_column = 0;
end

sino = (arg.phi_ob' * sino.').';
sino = arg.r_ob' * sino;
if arg.flag180
	sino = sino(:,1:end/2) + flipdim(sino(:,end/2+[1:end/2]),1);
end

if flag_column
	sino = reshape(sino, prod(arg.dimi), []);
end


%
% rebin_par2fan_test()
%
function rebin_par2fan_test()

gm = sino_geom('moj', 'nb', 28, 'na', 11, ...
	'dx', 1.1, 'offset_r', 0.2, 'orbit', 180);
gp = sino_geom('par', 'nb', 20, 'na', 11, ...
	'dr', 1, 'offset_r', 0.1, 'orbit', 180);
gf = sino_geom('fan', 'nb', 21, 'na', 9, ...
	'ds', 0.5, 'offset_s', 0.25, 'orbit', 360, ...
	'dsd', 949, 'dod', 408);

oo{1} = rebin_par2fan(zeros(gp.nb,gp.na), gp, gf, 'ob', 1);
oo{2} = rebin_par2fan(zeros(gm.nb,gm.na), gm, gf, 'ob', 1);
oo{3} = rebin_par2fan(zeros(gm.nb,gm.na), gm, gp, 'ob', 1);

for ii=1:3
	ob = oo{ii};

	% test adjoint
	A = ob(:,:);
	B = ob';
	B = B(:,:);
	if max_percent_diff(A, B') > 3e-5
		error 'adjoint'
	else
		printm 'adjoint ok'
	end
end
