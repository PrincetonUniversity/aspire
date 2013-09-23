  function ob = fatrix2_mtimes2(ob1, ob2, varargin)
%|function ob = fatrix2_mtimes2(ob1, ob2, options)
%|
%| Construct fatrix2 object that is the product of two objects:
%|	ob = ob1 * ob2.
%| Requires size(ob1,2) == size(ob2,1) as in matrix multiplication.
%|
%| in
%|	ob1	*atrix		any object that can do "mtimes" and "size"
%|	ob2	*atrix		''
%|
%| out
%|	ob	fatrix2		ob1 * ob2
%|
%| caution: this will have limited (if any) "block" capabilities
%| that are inherited from ob1 (if it has any):
%| 	ob{iblock} = ob1{iblock} * ob2		todo: under construction
%|
%| Copyright 2010-12-04, Jeff Fessler, University of Michigan

if size(ob1, 2) ~= size(ob2, 1)
	pr size(ob1,2)
	pr size(ob2,1)
	error 'size mismatch'
end
if isfield(ob1, 'idim') && isfield(ob2, 'odim') ...
	&& ~isequal(ob1.idim, ob2.odim)
	error 'idim/odim mismatch'
end

if isa(ob2, 'fatrix2')
	idim = ob2.idim;
	imask = ob2.imask;
else
	idim = size(ob2,2); 
	imask = [];
end

if isa(ob1, 'fatrix2') % && (numel(idim) > 1 || ~isempty(imask)) % caution
	odim = ob1.odim;
	omask = ob1.omask;
else
	odim = size(ob1,1);
	omask = [];
end

if isfield(ob2, 'nblock') && ~isempty(ob2.nblock)
	fail 'block object on rhs unsupported'
end

if isfield(ob1, 'nblock') && ~isempty(ob1.nblock)
	fail 'block object on lhs could be added - ask jf'
%	argblock = {'blockify_data', ob1.blockify_data, ...
%		'block_setup', @fatrix2_block_setup, ...
%		'mtimes_block', @fatrix2_mtimes_block};
end

if isa(ob1, 'fatrix2') && isa(ob2, 'fatrix2')
	forw = @(arg, x) (arg.ob1.scale * arg.ob2.scale) * ...
		arg.ob1.handle_forw(arg.ob1.arg, ...
		arg.ob2.handle_forw(arg.ob2.arg, x));
	back = @(arg, y) conj(arg.ob1.scale * arg.ob2.scale) * ...
		arg.ob2.handle_back(arg.ob2.arg, ...
		arg.ob1.handle_back(arg.ob1.arg, y));
else
	forw = @(arg, x) arg.ob1 * (arg.ob2 * x);
	back = @(arg, y) arg.ob2' * (arg.ob1' * y);
	if isa(ob1, 'fatrix2') || isa(ob2, 'fatrix2')
		warn('multiplying %s * %s may fail', class(ob1), class(ob2))
	end
end

arg.ob1 = ob1;
arg.ob2 = ob2;

ob = fatrix2('arg', arg, ...
	'idim', idim, 'imask', imask, ...
	'odim', odim, 'omask', omask, ...
	'forw', forw, 'back', back, ...
	'power', @fatrix2_mtimes2_power);
%	arg.block{:},


% fatrix2_mtimes2_power(): A.^p
function ob = fatrix2_mtimes2_power(ob, pow)
ob1 = ob.arg.ob1;
if isfield(ob1, 'caller') && streq(ob1.caller, 'diag', 4)
	ob.arg.ob1 = ob.arg.ob1 .^ pow;
	ob.arg.ob2 = ob.arg.ob2 .^ pow;
else
	error 'power defined only for diag * object'
end
