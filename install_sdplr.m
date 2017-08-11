% INSTALL_SDPLR Install the SDPLR package
%
% Usage
%	 install_sdplr();
%	 install_sdplr(url, location);
%
% Input
%	 url: The URL from which the package should be downloaded. By default,
%		this is
%
%		   http://sburer.github.io/files/SDPLR-1.03-beta.zip
%
%	 location: The location in which the package should be installed. By
%		default, this is the subfolder 'extern' of the ASPIRE root folder.
%
% Description
%	 This function downloads the SDPLR package by Burer, Monteiro, and Choi,
%	 which is a package for solving large semidefinite optimization problems.
%	 The package includes optimized code in C which must be compiled. For
%	 this reason, a working C compiler must be installed.

function install_sdplr(url, location)
	if nargin < 1 || isempty(url)
		url = 'http://sburer.github.io/files/SDPLR-1.03-beta.zip';
	end

	if nargin < 2 || isempty(location)
		location = fullfile(aspire_root(), 'extern');
	end

	fprintf('Installing the SDPLR package.\n');
	fprintf('URL: %s\n', url);
	fprintf('Location: %s\n', location);
	fprintf('\n');

	ind = find(url=='/', 1, 'last');
	filename = url(ind+1:end);

	filepath = fullfile(location, filename);

	if exist(filepath, 'file')
		fprintf('Package exists on disk. Skipping download.\n');
	else
		if ~exist(location, 'dir')
			mkdir(location);
		end

		try
			fprintf('Downloading...');
			urlwrite(url, filepath);
			fprintf('OK\n');
		catch
			fprintf('Failed.\n');
			fprintf('Please download the package at the URL\n');
			fprintf('	 %s\n', url);
			fprintf('store it in the location\n');
			fprintf('	 %s\n', filepath);
			fprintf('and rerun the installation.\n');

			return;
		end
	end

	% Save current dir so that we can go back.
	current_dir = pwd();

	% Move into folder since untar depends on current directory.
	cd(fileparts(filepath));

	unzipped_files = unzip(filepath);

	if numel(unzipped_files) == 0 || numel(unzipped_files{1}) == 0
		error(['No files extracted. The package is most likely already ' ...
			'installed.' char(10) 'To reinstall, please remove the ' ...
			'associated directory.']);
	end

	sdplr_dir = fileparts(unzipped_files{1});
	ind = find(sdplr_dir==filesep, 1);
	if ~isempty(ind)
		sdplr_dir = sdplr_dir(1:ind-1);
	end

	sdplr_root = fullfile(fileparts(filepath), sdplr_dir);

	% Remove the SDPLR root if it's in the path. Messing around with MEX files
	% that are potentially loaded can cause trouble.
	paths = strsplit(path, pathsep);
	if any(strcmp(sdplr_root, paths))
		rmpath(sdplr_root);
	end

	cd(sdplr_root);

	compile_sdplr();

	try
		x = sdplr([1 1], 1, [-2 -1]', struct('l', 2));
		fprintf('MEX files are working.\n');
	catch
		error('MEX compilation failed.');
	end

	cd(current_dir);

	addpath(sdplr_root);

	fprintf('SDPLR successfully installed!\n');
end

function compile_sdplr()
	% Adapted from 'mexinstall.m' of the SDPLR package with additions for
	% Octave and clang support.

	if isoctave()
		% For the source to compile under Octave, we need to make a few
		% changes.

		% Read in mexsdplr.c and split into lines.
		source_file = fileread(fullfile('source', 'mexsdplr.c'));
		source_lines = strsplit(source_file, '\n');

		% Add stddef.h to includes.
		ind = find(strcmp(source_lines, '#include <stdlib.h>'), 1);
		source_lines(ind+1:end+1) = source_lines(ind:end);
		source_lines{ind} = '#include <stddef.h>';

		% Remove matrix.h from includes. 
		ind = find(strcmp(source_lines, '#include "matrix.h"'), 1);
		source_lines(ind) = [];

		source_file = strjoin(source_lines, '\n');

		% Write modified file back to mexsdplr.c.
		fid = fopen(fullfile('source', 'mexsdplr.c'), 'w');
		fprintf(fid, '%s', source_file);
		fclose(fid);
	end

	% For things to compile under clang, we need to ensure that the main
	% function definition is correct.

	% Read in proto.h and split into lines.
	source_file = fileread(fullfile('source', 'proto.h'));
	source_lines = strsplit(source_file, '\n');

	% Change definition of main in proto.h.
	old_line = 'int main(size_t argc, char *argv[]);';
	new_line = 'int main(int argc, char *argv[]);';
	ind = find(strcmp(source_lines, old_line), 1);
	source_lines{ind} = new_line;

	source_file = strjoin(source_lines, '\n');

	% Write modified file back to proto.h.
	fid = fopen(fullfile('source', 'proto.h'), 'w');
	fprintf(fid, '%s', source_file);
	fclose(fid);

	% Read in main.c and split into lines.
	source_file = fileread(fullfile('source', 'main.c'));
	source_lines = strsplit(source_file, '\n');

	% Change definition of main in main.c.
	old_line = 'int main(size_t argc, char *argv[])';
	new_line = 'int main(int argc, char *argv[])';
	ind = find(strcmp(source_lines, old_line), 1);
	source_lines{ind} = new_line;

	source_file = strjoin(source_lines, '\n');

	% Write modified file back to main.c.
	fid = fopen(fullfile('source', 'main.c'), 'w');
	fprintf(fid, '%s', source_file);
	fclose(fid);

	winlib1 = [matlabroot '\extern\lib\win32\lcc\libmwlapack.lib'];
	winlib2 = [matlabroot '\extern\lib\win32\lcc\libmwblas.lib'];

	if ~isoctave()
		matlabversion = sscanf(version,'%f');
		matlabversion = matlabversion(1);
		matlabversion_1 = floor(matlabversion);
		matlabversion_2 = matlabversion - matlabversion_1;
		matlabversion_2 = sscanf(num2str(matlabversion_2), '0.%d\n');
	end

	subdir = 'source';

	fname = {'mexsdplr.c', ...
		'copystructures.c', ...
		'dataoper.c', ...
		'eigen.c', ...
		'initialize.c', ...
		'lbfgs.c', ...
		'linesearch.c', ...
		'main.c', ...
		'misc.c', ...
		'params.c', ...
		'rankreduce.c', ...
		'readdata.c', ...
		'sdplrlib.c', ...
		'timefuncs.c', ...
		'util.c'};

	if ~isoctave()
		mexcmd = 'mex -O CFLAGS="\$CFLAGS -std=iso9899:1999" -D__MEX ';
	else
		mexcmd = 'mex -O -std=iso9899:1999 -D__MEX ';
	end

	if ispc()
		mexcmd = [mexcmd '-D__WIN32 '];
	end

	if ~isoctave() && matlabversion_1 >= 7 && matlabversion_2 >= 3
		mexcmd = [mexcmd '-largeArrayDims '];
	elseif isoctave()
		% TODO: How to detect whether -largeArraydims is available in Octave?
	end

	for k = 1:numel(fname)
		mexcmd = [mexcmd fullfile(subdir, fname{k}) ' '];
	end

	mexcmd = [mexcmd fullfile('gsl-1.5', 'poly', 'eval.c') ' '];
	mexcmd = [mexcmd fullfile('gsl-1.5', 'poly', 'solve_cubic.c') ' '];

	if ispc() && ~isoctave()
		mexcmd = [mexcmd '"' WINLIB1 '" "' WINLIB2 '" '];
	elseif ~isoctave()
		mexcmd = [mexcmd '-lmwlapack -lmwblas '];
	elseif isoctave()
		mexcmd = [mexcmd strtrim(mkoctfile('-p', 'LAPACK_LIBS')) ' '];
		mexcmd = [mexcmd strtrim(mkoctfile('-p', 'BLAS_LIBS')) ' '];
	end

	eval(mexcmd);
end
