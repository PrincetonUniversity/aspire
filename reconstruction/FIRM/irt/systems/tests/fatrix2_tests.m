  function fatrix2_tests(A, varargin)
%|function fatrix2_tests(A, [option])
%|
%| perform basic tests of fatrix2 objects
%| in
%|	A	[nd np]		fatrix2 object
%|
%| option
%|	'x'	[idim]		optional input test image
%|	'name'	char		name of input object (optional for display)
%|	'caller' char		name of caller (optional for display)
%|	'full'	0|1		test A(:,:) and full(A)? (default: 1 if np<100)
%|	'halt'	0|1		if error, 1 (default) halt, 0 keyboard
%|	'complex' 0|1		if 1, test for complex inputs. default: 0
%|
%| Copyright 2005-8-2, Jeff Fessler, University of Michigan

if nargin < 1, help(mfilename), error(mfilename), end
if streq(A, 'test'), fatrix2_tests_test, return, end

arg.x = [];
arg.name = '';
arg.caller = '';
arg.complex = false;
arg.full = []; % see below
arg.halt = true; % halt if error
arg.test_block = true; % set to false to avoid recursion from within
arg = vararg_pair(arg, varargin);

if isempty(arg.full)
	arg.full = size(A,2) < 100; % hopefully small enough do manage A(:,:)
end

if isempty(arg.name)
	arg.name = inputname(1);
	arg.caller = caller_name;
end

o1d = numel(A.odim) == 1 & isempty(A.omask); % 1d output with full omask

fatrix2_test_sub(A) % subsref
fatrix2_test_dim(A, A.imask_array, A.omask_array, o1d)
fatrix2_test_fun(A, A.imask_array, A.omask_array)
arg.x = fatrix2_test_x(arg.x, A.imask_array, arg.complex);
fatrix2_test_basic(A, A.imask_array, A.omask_array, o1d, arg.x, ...
	arg.name, arg.caller, arg.complex, arg.full, arg.halt)
if arg.test_block
	fatrix2_test_block(A, arg.x, arg.complex, o1d) % [,] [;] etc.
end


% fatrix2_test_sub()
% test subsref ob.?
function fatrix2_test_sub(A)
A.arg;
A.caller;
A.imask;
A.omask;
A.idim;
A.odim;
A.meth;
A.docs;
A.scale;
A.size;
A.imask_array;
A.omask_array;


% fatrix2_test_dim()
% test dimensions
function fatrix2_test_dim(A, imask, omask, o1d)

fatrix2_check_dim1(A.imask_array, A.idim)
fatrix2_check_dim1(A.omask_array, A.odim)
nd = size(A,1);
np = size(A,2);
jf_equal(size(A), [nd np])

jf_equal(np, sum(imask(:)))
jf_equal(nd, sum(omask(:)))

x = single(imask);
tmp = A * x;
fatrix2_check_dim1(tmp, A.odim)

y = single(omask);
tmp = A' * y;
if o1d
	fatrix2_check_dim1(tmp, size(A,2))
else
	fatrix2_check_dim1(tmp, A.idim)
end

tmp = A * x(imask);
jf_equal(size(tmp), [nd 1])

tmp = A' * y(omask);
jf_equal(size(tmp), [np 1])

% check forw/back w.r.t. masking
tmp = A * ones([A.idim 1]);
if any(tmp(~omask))
	fail 'nonzero outside omask'
end

tmp = A' * ones([A.odim 1]);
if o1d
	jf_equal(size(tmp), [np 1])
else
	if any(tmp(~imask))
		fail 'nonzero outside imask'
	end
end


% fatrix2_test_fun()
function fatrix2_test_fun(A, imask, omask)

if 1 % multi-index: A(:,jj) A(ii,:)
	jj = [1 size(A,2)];
	tmp0 = A(:,jj);
	tmp0 = full(tmp0);
	tmp1 = A(:,jj(1));
	tmp2 = A(:,jj(2));
%	jf_equal(tmp0, [tmp1 tmp2])
	equivs(tmp0, [tmp1 tmp2])

	ii = [1 size(A,1)];
	tmp0 = A(ii,:);
	tmp0 = full(tmp0);
	tmp1 = A(ii(1),:);
	tmp2 = A(ii(2),:);
	jf_equal(tmp0, [tmp1; tmp2])
end

nd = size(A,1);
np = size(A,2);

% mtimes
a1 = A * ones(np,1);
r1 = ones(1,nd) * A;
equivs(sum(a1), sum(r1))

b1 = A' * ones(nd,1);
c1 = ones(1,np) * A';
equivs(sum(b1), sum(c1))

if np < 100 % small enough for 'full'
	tmp = A(:,:);
	full(A);
	full(A');
%	jf_equal(A(end), tmp(end))
	equivs(A(end), tmp(end))
end

if 0
	B = A';
	B(:,2);
	C = A(:,[3 1]);
	jf_equal(C(:,1), A(:,3))
	C = B(:,[2 7]);
	C(:,2);
	C(:,1:2);
	full(C(:,1:2));
end


% fatrix2_test_x()
function x = fatrix2_test_x(x, imask, do_complex)

if do_complex && ~isempty(x) && isreal(x)
	fail 'asked for complex test but provided real x'
end

if ~isempty(x), return, end % defer to user's x

switch ndims(imask)
case 2
	[nx ny] = size(imask);
	ig = image_geom('nx', nx, 'ny', ny, 'dx', 1);
	x = ellipse_im(ig, 'shepplogan-emis', 'type', 'slow');
	clear ig nx ny

case 3
	[nx ny nz] = size(imask);
	ig = image_geom('nx', nx, 'ny', ny, 'nz', nz, 'dx', 1, 'dz', 1);
	x = ellipsoid_im(ig, []);
	clear ig nx ny nz

otherwise
	error 'only 2d and 3d implemented'
end

if do_complex
	x = x + 1i * flipdim(x, 1); % something complex for testing
end
x = x .* imask; % test with mask-limited object


% fatrix2_test_basic()
function fatrix2_test_basic(A, imask, omask, o1d, x, name, caller, ...
	do_complex, do_full, do_halt)

% sum
sum(sum(A));

% trick:
% some tests are meaningful only for objects that can go forward
% into an array instead of a vector, i.e., for tomography objects
% that make sinograms, but not for Gmri.

ya = A * x; % array
yv = A * x(imask); % vector

if ~isreal(ya) && ~do_complex
	warn 'recommend testing complex objects with ''complex'' option'
end

% scale
B = 7 * A;
y2 = B * x;
jf_equal(7 * ya, y2)

% A * : array vs vector
y2 = ya(omask);
if do_halt
	jf_equal(yv, y2)
else
	try
		jf_equal(yv, y2)
	catch
		warn 'A * x vs reshape(A * x(imask)) failed'
		keyboard
	end
end

% A * [x() x()]
if 1
	y2 = A * [x(imask) 2*x(imask)];
	z2 = [yv 2*yv];
	if do_halt
		equivs(y2, z2) % may not be *exactly* equal due to 2*
	else
		try
			equivs(y2, z2)
		catch
			warn 'A * [x() x()] failed'
			keyboard
		end
	end
end

% version of ndims that returns "1" for column vectors!
ndims1 = @(x) ndims(x) - (size(x,ndims(x)) == 1);
catnext = @(a,b) cat(1+ndims1(a),a,b);

%if isvar('A.arg.odim') % caution: handle single view projectors [ns nt 1 nrep]
%	catp = @(a,b) cat(1+length(A.arg.odim), a, b);
%else
	catp = catnext;
%end

% A * [x x]
if 1
	y2 = A * catnext(x, 2*x);
	z2 = catp(ya, 2*ya);
	if do_halt
		equivs(y2, z2) % may not be *exactly* equal due to 2*
	else
		try
			equivs(y2, z2)
		catch
			warn 'A * [x x] failed'
			keyboard
		end
	end
end

y0 = ya;
if do_complex && isreal(y0)
	warn 'complex testing is incomplete; tell jeff!'
end

% A' * : array vs vector
xa = A' * ya; % array
xv = A' * yv; % vector

if o1d
	xa = embed(xa, imask);
end

x2 = xa(imask);
jf_equal(xv, x2)

x2 = embed(xv, imask);
jf_equal(xa, x2)

% A' * [y() y()]
if 1
	x2 = A' * [yv 2*yv];
	v2 = [xv 2*xv];
	equivs(x2, v2) % may not be *exactly* equal due to 2*
end

% A' * [y y]
if 1
	tmp = catp(ya, 2*ya);
	x2 = A' * tmp;
	if o1d
		x2 = embed(x2, imask);
	end
	v2 = catnext(xa, 2*xa);
	equivs(x2, v2) % may not be *exactly* equal due to 2*
end

% indexing single and multiple arguments and test consistency thereof
% A(:,?)
c1 = A(:,1);
c2 = A(:,end);
c12 = A(:, [1 end]);
c12 = [c12(:,1) c12(:,2)]; % for fatrix2
jf_equal(c12, [c1 c2])

% A(?,:)
r1 = A(1,:);
r2 = A(end,:);
r12 = A([1 end], :);
r12 = [r12(1,:); r12(2,:)]; % for fatrix2
jf_equal(r12, [r1; r2])

printm('passed for "%s" in %s', name, caller)


if do_full
	Af = A(:,:);
	jj = 1;
	tmp = A(:,jj);
	equivs(tmp, Af(:,jj))
	ii = 1;
	tmp = A(1,:);
	jf_equal(tmp, Af(1,:))
end


% fatrix2_test_block()
% A1+A2
% [A1; A2]
% [A1, A2]
function fatrix2_test_block(A, x, do_complex, o1d)

test_fun = @(B) fatrix2_tests(B, 'complex', do_complex, 'x', x, ...
			'test_block', false); % trick: avoid recursion

y1 = A * x;
B = A + 7 * A;
y2 = B * x;
equivs(y2, 8*y1)
test_fun(B)

%B = [A, 7 * A, -A]; % todo horzcat

B = [A; 7 * A; -A]; % vertcat, same odim
y2 = B * x;
dimcat = B.arg.dim_cat;
yt = cat(dimcat, y1, 7*y1, -y1); 
try
	jf_equal(y2, yt)
catch
	keyboard
end
test_fun(B)

B = [A; [7 * A; -A]]; % vertcat, different odim
y2 = B * x;
yt = [y1(:); 7*y1(:); -y1(:)]; % trick: 1d due to different odims
try
	jf_equal(y2, yt)
catch
	keyboard
end
test_fun(B)

B = build_gram(A);
x1 = A' * (A * x); % array
if o1d
	x1 = embed(x1, A.imask);
end
x2 = B * x;
equivs(x2, x1)

y1 = A' * (A * x(A.imask_array)); % vector
y2 = B * x(A.imask_array);
equivs(y2, y1)
y2 = B' * x(A.imask_array);
equivs(y2, y1)
test_fun(B)


% fatrix2_mtimes2_test
% A * B
function fatrix2_mtimes2_test(A, B, x, do_complex)
C = A * (7 * B);
y1 = 7 * (A * (B * x));
y2 = C * x;
jf_equal(y1, y2)
fatrix2_tests(C, 'x', x, 'complex', do_complex)


% fatrix2_tests_test()
function fatrix2_tests_test

psf = [0 1 2 1 0; 1 2 4 3 1; 0 2 3 1 0];
psf = psf / sum(psf(:));
idim = [24 30];
mask = true(idim);
mask(1) = 0;
A = Gblur(mask, 'psf', psf);

x = fatrix2_test_x([], mask, false);

if 1
	fatrix2_tests(A, 'x', x)

	A = Gblur(mask, 'psf', psf + 1i);
	x = fatrix2_test_x([], mask, true);
	fatrix2_tests(A, 'x', x, 'complex', 1)

	fatrix2_mtimes2_test(A', A, x, true)
	fatrix2_mtimes2_test(A, A', x, true)
end

if 1
	B = Gdiag(ones(sum(mask(:)),1), 'mask', mask(:));
	T = A' * A;
	C = B + T; % stress test mismatched idim,odim
	y1 = B * x(mask) + A' * A * x(mask);
	y2 = C * x(mask);
	equivs(y1,y2)
	y3 = C' * x(mask);
	equivs(y1,y3)
end
