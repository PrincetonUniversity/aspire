  function ob = Gtomosyn(cg, ig, varargin)
%|function ob = Gtomosyn(cg, ig, [options])
%|
%| Construct Fatrix object for forward and back-projection with
%| a 3D tomosynthesis geometry having a flat detector and source moving in arc.
%| Currently is simply a "stack" of Gcone objects - not efficient for multicore.
%|
%| in
%|	cg	ct_geom object (with orbit=0, for head-on projection view)
%|	ig	image_geom object (recon grid dimensions, FOV and mask)
%|
%| option
%|	'angles'	list of projection angles in degrees. default: -30:2:30
%|	'type'		see Gcone.m
%|	'is_ns_nt'	1 (default) for [ns nt na] or 0 for [nt ns na]
%|	'nthread'	# of threads, default: jf('ncore')
%|	'scale'		scale mex output by this. (default: [])
%|
%| out
%|	ob [nd np]	np = sum(mask(:)), so it is already "masked"
%|			nd = [ns * nt * na] if is_ns_nt = 1
%|			nd = [nt * ns * na] if is_ns_nt = 0
%|
%| Typically for a mx * my detector, ns=nx=mx and nt=nz=my
%| For more help, see Gcone.m and Gtomosyn_test.m
%|
%| Copyright 2010-07-27, Jeff Fessler, University of Michigan

if nargin == 1 & streq(cg, 'test'), run_mfile_local Gtomosyn_test, return, end
if nargin < 2, help(mfilename), error(mfilename), end

% option defaults
arg.angles = -30:3:30;
arg.type = 'sf2';
arg.chat = false;
arg.scale = [];
arg.is_ns_nt = true; % default: [ns nt na] (channel fastest, row, view slowest)
arg.nthread = jf('ncore');

% options specified by name/value pairs
arg = vararg_pair(arg, varargin);

thetas = deg2rad(arg.angles);

if cg.orbit ~= 0 ||  cg.orbit_start ~= 0
	fail 'orbit'
end

na = length(thetas);
Ac = cell(na,1);
for ia=1:na
	theta = thetas(ia);
	[cg_tmp ig_tmp] = Gtomosyn_setup(cg, ig, theta);
	Ac{ia} = Gcone(cg_tmp, ig_tmp, ...
		'type', arg.type, ...
		'scale', arg.scale, ...
		'is_ns_nt', arg.is_ns_nt, ...
		'nthread', arg.nthread, ...
		'chat', arg.chat);
	if arg.chat
		cg_tmp.plot(ig_tmp);
		titlef('ia %d', ia)
		xlabel 'x', ylabel 'z !!'
	prompt
	end
end

ob = block_fatrix(Ac, 'type', 'col', 'tomo', true);


% Gtomosyn_setup()
function [cg_out ig_out] = Gtomosyn_setup(cg, ig, theta)
dsd = cg.dsd;
cg_out = cg;
ig_out = ig;

shift = dsd * sin(theta);

cg_out.dsd = dsd * cos(theta);
cg_out.offset_s = cg.offset_s + shift / cg.ds;
ig_out.offset_x = ig.offset_x + shift / ig.dx;
