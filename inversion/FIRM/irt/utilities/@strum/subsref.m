 function out = subsref(ob, args)
%function out = subsref(ob, args)
% handle subscript references like ob.ref or ob.ref(arg[s])
% Copyright 2006-1-19, Jeff Fessler, The University of Michigan

st = struct(ob);

if args(1).type ~= '.'
	error 'only ob.data and ob.meth(args) done'
end

%
% ob.name
%
name = args(1).subs;
if isfield(st, name)
	if length(args) == 1
		out = st.(name);
	else
		out = subsref(st.(name), args(2:end));
	end

%
% ob.name, where name is in data field
%
elseif isfield(st.data, name)
	out = st.data.(name);
	if length(args) > 1

		% ob.name{:} does not work properly because there seems
		% to be no way to return an arbitrary number of arguments
		% (in a list) to the caller.  this is unfortunate because
		% a.b = {1,2} and a.b{:} works for ordinary structures.
% todo: see this solution involving numel()
% http://www.mathworks.com/support/solutions/data/1-1ABOD.html?solution=1-1ABOD
		if iscell(out) && ...
		(isequal(args(2).subs, {':'}) || length(args(2).subs{1}) > 1)
			warn 'instead of strum.data{:} etc., use two steps:'
			warn 'tmp = strum.data; tmp{:} or tmp{1:end} etc.'
			fail 'strum.data{:} etc. unsupported thanks to matlab'
		end

		% fix: strum1.strum2.plot() does not work here
		out = subsref(out, args(2:end));
	end

%
% ob.name or ob.name(arg[s]), where name is in meth field
%
elseif isfield(st.meth, name)
	fun = st.meth.(name); % function handle for method
	if length(args) == 1
		if isfreemat || nargout(fun) % freemat: does not like nargout()
%			out = feval(fun, ob);
			out = fun(ob);
		else
%			feval(fun, ob);
			fun(ob);
		end
	elseif length(args) == 2
		if isfreemat || nargout(fun) % freemat: does not like nargout()
%			out = feval(fun, ob, args(2).subs{:});
			out = fun(ob, args(2).subs{:});
		else
%			feval(fun, ob, args(2).subs{:});
			fun(ob, args(2).subs{:});
		end
	else
		printm('unknown subs')
		keyboard
	end

else
	disp(st)
	fail(['unknown field name: ' name])
end
