  function x = read_ncat(varargin)
%|function x = read_ncat(options)
%|
%| read in a slice of the NCAT phantom
%| and assign it attenuation coefficients in inverse mm units.
%|
%| options
%|	'nx'		desired size
%|	'ny'
%|	'mu'	[5 1]	desired image intensity values
%|			0 background, 1 lung, 2 body, 3 spine, 4 ribs
%|	'marrow' 0|1	assign (fake) marrow to class 5 (default: 0)
%|	'ddir'		data directory
%| out
%|	x	[nx ny]	ncat image
%|
%| Jeff Fessler, University of Michigan

if nargin == 1 && streq(varargin{1}, 'test'), read_ncat_test, return, end
if ~nargout, help(mfilename), error(mfilename), end

arg.nx = 1024;
arg.ny = [];
arg.mu = [];
arg.marrow = false;
arg.marrow_erode = 6; % erode 6 bone pixels to make marrow
arg.ddir = '';
arg.file = 'ncat,1024,slice500,0to4.fld';
arg = vararg_pair(arg, varargin);
if isempty(arg.ny), arg.ny = arg.nx; end

x = read_ncat_do(arg.ddir, arg.file, arg.nx, arg.ny, arg.mu, ...
	arg.marrow, arg.marrow_erode);


%
% read_ncat_do()
%
function x = read_ncat_do(ddir, file, nx, ny, mu, do_marrow, marrow_erode)

% guess .../data directory by looking parallel to '.../transmission' directory
if ~isvar('ddir') | isempty(ddir)
	t = path_find_dir([filesep 'transmission']);
	ddir = strrep(t, 'transmission', 'data');
end

if ~exist(ddir, 'dir')
	warn('cannot find data directory %s', ddir)
	fail('edit path in %s.m for your installation', mfilename)
end

file = [ddir filesep 'ncat,1024,slice500,0to4.fld'];
x = fld_read(file);

if do_marrow
	tmp = x >= 3; % spine and ribs
	marrow = imerode(tmp, strel('disk', marrow_erode));
	x(marrow) = 5; % new class
end

if ~isempty(mu)
	x = mu(1+x); % map indices to desired image intensity values
end

x = phantom_resize(x, nx, ny);


function read_ncat_test
% original way:
%x = read_ncat('nx', 512, 'ny', 450, 'mu', [0 0.05 0.2 0.3 0.4]/0.2*1000);
% new way with marrow:
x = read_ncat('nx', 512, 'ny', 450, 'marrow', true, ...
	'mu', [0 0.05 0.2 0.4 0.4 0.2]/0.2*1000);
im(x, 'ncat phantom'), cbar
