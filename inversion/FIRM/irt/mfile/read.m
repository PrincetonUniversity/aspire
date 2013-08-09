 function test = read(arg1, arg2, arg3, arg4, arg5)
%function test = read(dir, file, slice, frame, chat)
%function test = read(file, slice, frame, chat)
% read a file in directory 'dir' (which is optional)
% or optionally one slice,frame of a file

slice = [];
frame = 0;
chat = 0;

if nargin < 1 | nargin > 5
	help read
	error args
elseif nargin == 1
	file = arg1;
else
	if isstr(arg2)
		file = [arg1 arg2];
		if nargin >= 3,	slice	= arg3;	end
		if nargin >= 4,	frame	= arg4;	end
		if nargin >= 5,	chat	= arg5;	end
	else
		file = arg1;
		if nargin >= 2,	slice	= arg2;	end
		if nargin >= 3,	frame	= arg3;	end
		if nargin >= 4,	chat	= arg4;	end
	end
end

if length(file) > 4
	ismat = strcmp(file((end-3):end), '.mat');
else
	ismat = 0;
end

	op = sprintf('op -chat %d ', chat);

	if exist(file) ~= 2
%		error(['File ' file ' not found!'])
	end

	tmp = '/tmp/tmp.mat';
	if exist(tmp) == 2
		eval(['!/bin/rm -f ' tmp])
	end

%	type = 'double';
	if ismat
		type = 'double';
	else
%		type = '-';
		type = 'raw';	% ask op for raw file format
	end

	if ~isempty(slice)
		if length(slice) == 1
			choose = sprintf('%d %d %d %d 0 0', slice, slice, frame, frame);
		elseif length(slice) == 2
			choose = sprintf('%d %d %d %d 0 0', slice(1), slice(2), frame, frame);
		else
			error slice
		end
		comm = [op 'slice ' tmp ' ' file ' ' type ' ' choose];
	else
		comm = [op 'convert ' tmp ' ' file ' ' type ' test'];
	end

	if chat
		disp(comm)
	end
	os_run(comm)
	if ~exist(tmp, 'file')
		error(['in read.m: Conversion error on: ' file])
	end

	load(tmp)
	eval(['!/bin/rm -f ' tmp])
