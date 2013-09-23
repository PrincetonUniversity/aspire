  function ob = Gblur(mask, varargin)
%|function ob = Gblur(mask, options)
%|
%| Construct Gblur object for image restoration.
%|
%| See Gblur_test.m for example usage.
%|
%| in
%|	mask	size(image)	logical array of object support.
%|
%| options
%|	'chat'		verbose printing of debug messages
%|	'psf'		point spread function (aka impulse response)
%|	'type'		type of blur:
%|				'conv,same'	usual case (default)
%|				'fft,same'	(periodic end conditions)
%|				todo: allow replicated end conditions!
%|				'imfilter,same'	(requires image toolbox)
%|				'imfilter,circ'	(requires image toolbox)
%|	'imfilter_options'	options to 'imfilter' (if used)
%|
%| out
%|	ob [nd np]	np = sum(mask(:)), so it is already "masked"
%|			nd = prod(size(mask)) for 'conv,same' type
%|
%| Copyright 2005-4-22, Jeff Fessler, University of Michigan

if nargin == 1 && streq(mask, 'test'), Gblur_test, return, end
if nargin < 1, help(mfilename), error(mfilename), end

arg.mask = mask;

% option defaults
arg.chat = 0;
arg.psf = 1; % identity matrix
arg.type = 'conv,same';
arg.class = 'fatrix2';
arg.imfilter_options = {};

% options specified by name/value pairs
arg = vararg_pair(arg, varargin);

if ndims(arg.psf) ~= ndims(mask), error 'psf dim mismatch', end
if any(size(arg.psf) > size(mask)) error 'psf too large', end

if streq(arg.type, 'imfilter,circ') && exist('imfilter') ~= 2
	warn 'no imfilter (no image processing toolbox) so using fft'
	arg.type = 'fft,same';
end

if streq(arg.type, 'imfilter,same') && exist('imfilter') ~= 2
	warn 'no imfilter (no image processing toolbox) so using slower convn'
	arg.type = 'conv,same';
end

does_many = false;
switch arg.type
case 'conv,same'
	arg.odim = size(mask);
	arg.psf_flip = conj(flipdims(arg.psf));
	handle_forw = @(arg,x) convn(x, arg.psf, 'same');
	handle_back = @(arg,y) fatrix2_maskit(arg.mask, ...
			convn(y, arg.psf_flip, 'same'));
	does_many = true;

case 'fft,same'
	arg.odim = size(mask);
	[arg handle_forw handle_back] = Gblur_setup_fft_same(arg);

case 'imfilter,circ'
	arg.odim = size(mask);
	arg.psf_flip = conj(flipdims(arg.psf));
	handle_forw = @(arg,x) imfilter(x, arg.psf, 'conv', 'same', ...
		'circular', arg.imfilter_options{:});
	handle_back = @(arg,y) fatrix2_maskit(arg.mask, ...
			imfilter(y, arg.psf_flip, 'conv', 'same', ...
				'circular', arg.imfilter_options{:}));

case 'imfilter,same'
	arg.odim = size(mask);
	arg.psf_flip = conj(flipdims(arg.psf));
	handle_forw = @(arg,x) imfilter(x, arg.psf, 'conv', 'same', ...
		arg.imfilter_options{:});
	handle_back = @(arg,y) fatrix2_maskit(arg.mask, ...
			imfilter(y, arg.psf_flip, 'conv', 'same', ...
				arg.imfilter_options{:}));

otherwise
	error 'unknown blur type'
end

switch arg.class

case 'Fatrix'

	arg.nd = prod(arg.odim);
	arg.np = sum(mask(:));
	dim = [arg.nd arg.np]; % trick: make it masked by default!

	ob = Fatrix(dim, arg, 'caller', 'Gblur', ...
		'forw', @Gblur_forw_Fatrix, 'back', @Gblur_back_Fatrix);

case 'fatrix2'

	ob = fatrix2('idim', arg.odim, 'arg', arg, 'does_many', does_many, ...
		'forw', handle_forw, 'back', handle_back);

otherwise
	fail('bug')
end


% Gblur_setup_fft_same()
function [arg handle_forw handle_back] = Gblur_setup_fft_same(arg)

mask = arg.mask;

% put psf in center of array
tmp = zeros(size(mask));
if any(~rem(size(arg.psf),2)), error 'psf size must be odd', end
switch ndims(arg.psf)
case 2
	h1 = (size(arg.psf,1)-1)/2;
	h2 = (size(arg.psf,2)-1)/2;
	i1 = floor(size(mask,1)/2) + 1 + [-h1:h1];
	i2 = floor(size(mask,2)/2) + 1 + [-h2:h2];
	tmp(i1,i2) = arg.psf;
case 3
	h1 = (size(arg.psf,1)-1)/2;
	h2 = (size(arg.psf,2)-1)/2;
	h3 = (size(arg.psf,3)-1)/2;
	i1 = floor(size(mask,1)/2) + 1 + [-h1:h1];
	i2 = floor(size(mask,2)/2) + 1 + [-h2:h2];
	i3 = floor(size(mask,3)/2) + 1 + [-h3:h3];
	tmp(i1,i2,i3) = arg.psf;
otherwise
	error 'only 2d & 3d done'
end
arg.psf_fft = fftn(fftshift(tmp));

handle_forw = @(arg,x) ifftn(fftn(x) .* arg.psf_fft);
handle_back = @(arg,y) fatrix2_maskit(arg.mask, ...
			ifftn(fftn(y) .* conj(arg.psf_fft)));


% Gblur_forw_Fatrix(): y = A * x
function y = Gblur_forw_Fatrix(arg, x)

[x ei] = embed_in(x, arg.mask, arg.np);

switch arg.type
case 'conv,same'
	y = convn(x, arg.psf, 'same');
case 'fft,same'
	if ndims(x) > ndims(arg.mask);
		y = fft_conv_multi(x, arg.psf_fft);
	else
		y = ifftn(fftn(x) .* arg.psf_fft);
	end
case 'imfilter,circ'
	y = imfilter(x, arg.psf, 'conv', 'same', ...
		'circular', arg.imfilter_options{:});
case 'imfilter,same'
	y = imfilter(x, arg.psf, 'conv', 'same', arg.imfilter_options{:});
otherwise
	error 'bug'
end

y = ei.shape(y);


% Gblur_back_Fatrix(): x = A' * y
function x = Gblur_back_Fatrix(arg, y)

[y eo] = embed_out(y, arg.odim);

switch arg.type
case 'conv,same'
	x = convn(y, arg.psf_flip, 'same');
case 'fft,same'
	if ndims(y) > ndims(arg.mask);
		x = fft_conv_multi(y, conj(arg.psf_fft));
	else
		x = ifftn(fftn(y) .* conj(arg.psf_fft));
	end
case 'imfilter,circ'
	x = imfilter(y, arg.psf_flip, 'conv', 'same', ...
		'circular', arg.imfilter_options{:});
case 'imfilter,same'
	x = imfilter(y, arg.psf_flip, 'conv', 'same', arg.imfilter_options{:});
otherwise
	error 'bug'
end

x = eo.shape(x, arg.mask, arg.np);


%
% fft_conv_multi()
%
function y = fft_conv_multi(x, H);
dim = size(x);
y = zeros(size(x));
y = [];
for ii=1:dim(end)
	tmp = ifftn(fftn(stackpick(x, ii)) .* H);
	y = cat(length(dim), y, tmp); % slow: keeps enlarging y.
end
