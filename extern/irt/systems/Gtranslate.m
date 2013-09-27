  function ob = Gtranslate(mask, varargin)
%|function ob = Gtranslate(mask, options)
%|
%| Construct Gtranslate object for image registration.
%| (Internally stores image-sized arrays so not very memory efficient.)
%|
%| See Gtranslate_test() below for example usage.
%|
%| in
%|	mask	size(image)	logical array of object support.
%|
%| options
%|	'shift'	[ndim 1]	shift amount (default 0)
%|	'type'	char		translation method:
%|				'bspline3' todo (adjoint not working)
%|				'circshift' (integer shifts only)
%|				'fft' (terrible for non-integer shifts)
%|				'interpn,linear' linear interpolation (default)
%|
%| out
%|	ob [nd np]	np = sum(mask(:)), so it is already "masked"
%|			nd = np for 'conv,same' type
%|
%| Copyright 2010-3-22, Jeff Fessler, University of Michigan

if nargin == 1 && streq(mask, 'test'), Gtranslate_test, return, end
if nargin < 1, help(mfilename), error(mfilename), end

arg.mask = mask;

% option defaults
arg.shift = 0; % identity matrix
arg.type = 'interpn,linear';

% options specified by name/value pairs
arg = vararg_pair(arg, varargin);

arg.ndim = ndims(mask);
if arg.ndim == 2 && size(mask,2) == 1
	arg.ndim = 1;
end
if length(arg.shift) ~= arg.ndim
	fail 'ndim mismatch'
end
arg.np = sum(mask(:));
arg.nd = prod(size(mask));
dim = [arg.nd arg.np];

switch arg.type

case 'bspline3'
	if ndims(mask) ~= 2, fail 'only 2d done', end

	ig = image_geom('nx', size(mask,1), 'dx', 1, ...
			'ny', size(mask,2), 'dy', 1);
	kg = knot_geom( 'nx', 1, 'mx', 8*ig.nx, 'offset_x', ig.nx/2, ...
			'ny', 1, 'my', 8*ig.ny, 'offset_y', ig.ny/2);
	Bx = makeB(ig, kg);
	By = Bx;
	alphax = zeros(kg.dim, 'single');
	alphay = zeros(kg.dim, 'single');
	alphax(1) = -2.25 * arg.shift(1); % trick: empirical
	alphay(1) = -2.25 * arg.shift(2);
	W = makeW({Bx, By}, {alphax, alphay});

	arg.fun_forw = @(x) W * BsplVal2CoMirr(single(x));
	arg.fun_back = @(y) W' * y; % todo:
%	arg.fun_back = @(y) BsplCo2ValTranMirr(W' * y);

case 'circshift'
	if any(round(arg.shift) ~= arg.shift) % need to use fft?
		fail 'circshift needs integer shifts'
	end
	arg.fun_forw = @(x) circshift(x, arg.shift);
	arg.fun_back = @(y) arg.mask .* circshift(y, -arg.shift);

case 'fft'
	Nd = size(mask);
	for id = 1:arg.ndim
		N = Nd(id);
		kk = single([0:N-1]);
		phase{id} = exp(-2i * pi / N * kk * arg.shift(id));
	end
	tmp = ndgrid_jf('mat', phase);
	tmp = prod(tmp, 1 + arg.ndim);
	arg.phase = single(tmp);

%	arg.fun_forw = @(x) ifftn(fftn(x) .* phase);
	arg.fun_forw = @(x) translate_fft(x, arg.phase);
	arg.fun_back = @(y) arg.mask .* translate_fft(y, conj(arg.phase));

case 'interpn,linear'
	Nd = size(mask);
	for id = 1:arg.ndim
		N = Nd(id);
		xf{id} = single([1:N]) - arg.shift(id);
		xb{id} = single([1:N]) + arg.shift(id);
	end
	xf = ndgrid_jf('cell', xf);
	xb = ndgrid_jf('cell', xb);
	arg.fun_forw = @(x) interpn(x, xf{:}, 'linear', 0);
	arg.fun_back = @(x) interpn(x, xb{:}, 'linear', 0);

otherwise
	fail('type "%s" unknown', arg.type)
end

%
% build Fatrix object
%
ob = Fatrix(dim, arg, 'caller', 'Gtranslate', ...
	'forw', @Gtranslate_forw, 'back', @Gtranslate_back);


%
% Gtranslate_forw(): y = A * x
%
function y = Gtranslate_forw(arg, x)

[x ei] = embed_in(x, arg.mask, arg.np);

if ndims(x) > arg.ndim
	y = zeros(arg.nd, 1, class(x));
	for ll=1:size(x, ndims(x))
		y(:,ll) = col(arg.fun_forw(stackpick(x, ll)));
	end
	y = reshapee(y, size(arg.mask), []);
else
	y = arg.fun_forw(x);
end

y = ei.shape(y);


%
% Gtranslate_back(): x = A' * y
% (adjoint)
%
function x = Gtranslate_back(arg, y)

[y eo] = embed_out(y, size(arg.mask));

if ndims(y) > ndims(arg.mask)
	x = zeros(arg.nd, 1, class(y)); % trick!
	for ll=1:size(y, ndims(y))
		x(:,ll) = col(arg.fun_back(stackpick(y, ll)));
	end
	x = reshapee(x, size(arg.mask), []);
else
	x = arg.fun_back(y);
end

x = eo.shape(x, arg.mask, arg.np);


%
% translate_fft()
%
function y = translate_fft(x, phase);
y = ifftn(fftn(x) .* phase);
if isreal(x)
	y = real(y);
end


%
% Gtranslate_test()
%
function Gtranslate_test
mask = true(8,1); ishift = 2; fshift = [2.3]; % 1d
mask = true(10,8); ishift = [7 3]; fshift = [2.3 1.7]; % 2d
mask(1) = false;
A = Gtranslate(mask, 'type', 'circshift', 'shift', ishift);
x = mask .* rand(size(mask));
jf_equal(x, A' * (A * x))
jf_equal(x(mask), A' * (A * x(mask)))
Fatrix_test_basic(A, mask)

A = Gtranslate(mask, 'type', 'fft', 'shift', fshift);
Fatrix_test_basic(A, mask, 'complex', 0)
Fatrix_test_basic(A, mask, 'complex', 1)

A = Gtranslate(mask, 'type', 'interpn,linear', 'shift', fshift);
Fatrix_test_basic(A, mask)
Fatrix_test_basic(A, mask, 'complex', 1)

%A = Gtranslate(mask, 'type', 'bspline3', 'shift', fshift);
%Fatrix_test_basic(A, mask) % not working

if 0
	tmp = 0 * mask; tmp(ceil(end/2), ceil(end/2)) = 1;
	im plc 2 3
	im(1, tmp)
	tmp = A * tmp;
	im(2, abs(tmp), 'abs'), cbar
	im(3, real(tmp), 'real'), cbar
	im(4, imag(tmp), 'imag'), cbar
return
end

% try it for image registration
%f.type = 'fft'; % bad results!
%f.type = 'bspline3';
f.type = 'interpn,linear';

ig = image_geom('nx', 128, 'dx', 1);
f1 = ellipse_im(ig, 'shepplogan-emis', 'oversample', 2);

im plc 1 2
im(1, f1), cbar

mask = true(size(f1));
A = Gtranslate(mask, 'type', f.type, 'shift', fshift);
y1 = A * f1;


if 0 % examine integer shifts
	Ac = Gtranslate(mask, 'type', 'circshift', 'shift', ishift);
	yc = Ac * f1;
	Ab = Gtranslate(mask, 'type', f.type, 'shift', ishift);
	yb = Ab * f1;
	yb = reshape(yb, size(mask));
	im clf, im_toggle(yb, yc, [0 255])
return
end


if 0 % examine shifts
	Af = Gtranslate(mask, 'type', f.type, 'shift', fshift);
	yf = Af * f1;
	Ab = Gtranslate(mask, 'type', 'bspline3', 'shift', fshift);
	yb = Ab * f1;
	yb = reshape(yb, size(mask));
	im clf, im_toggle(yf, yb, [0 255])
return
end


shifts = linspace(0, 4, 21);
cost = zeros(size(shifts));
for is=1:length(shifts)
	tmp = [shifts(is) fshift(2)]; % in dim1 only
	A2 = Gtranslate(mask, 'type', f.type, 'shift', tmp);
	y2 = A2 * f1;
	cost(is) = norm(y2(:) - y1(:))^2;
end

if im
	im subplot 2
	plot(shifts, cost, '-o', fshift(1), 0, 'rx')
end
