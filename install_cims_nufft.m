% INSTALL_CIMS_NUFFT Install CIMS NUFFT library
%
% Usage
%    install_cims_nufft();
%    install_cims_nufft(url, location, force_compile);
%
% Input
%    url: The url from which the package should be downloaded. By default,
%       this is
%
%          http://cims.nyu.edu/cmcl/nufft/nufftall-1.3.3.tar.gz
%
%    location: The location in which the package should be installed. By
%       default, this is the subfolder 'extern' of the ASPIRE root folder.
%    force_compile: If set to yes, the precompiled MEX files shipped with the
%       package are removed and recompiled (default false).
%
% Description
%    This function downloads the NUFFT library from CIMS by Greengard & Lee,
%    which implements non-uniform fast Fourier transforms. Since these depend
%    on optimized Fortran code, the function checks whether the precompiled
%    versions that ship with the package work, and if not, it recompiles them.
%    For this task, the 'gfortran' compiler must be installed.

function install_cims_nufft(url, location, force_compile)
	if nargin < 1 || isempty(url)
		url = 'http://cims.nyu.edu/cmcl/nufft/nufftall-1.3.3.tar.gz';
	end

	if nargin < 2 || isempty(location)
		location = fullfile(aspire_root(), 'extern');
	end

	if nargin < 3 || isempty(force_compile)
		force_compile = false;
	end

	fprintf('Installing the CIMS NUFFT package.\n');
	fprintf('URL: %s\n', url);
	fprintf('Location: %s\n', location);
	fprintf('Force recompilation: %d\n', force_compile);
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
			fprintf('    %s\n', url);
			fprintf('store it in the location\n');
			fprintf('    %s\n', filepath);
			fprintf('and rerun the installation.\n');

			return;
		end
	end

	% Save current dir so that we can go back.
	current_dir = pwd();

	% Move into folder since untar depends on current directory.
	cd(fileparts(filepath));

	unzipped_files = gunzip(filepath);
	if numel(unzipped_files) == 1 && strcmp(unzipped_files{1}(end-3:end), '.tar')
		untarred_files = untar(unzipped_files{1});
		delete(unzipped_files{1});
		unzipped_files = untarred_files;
	end

	nufft_root = fullfile(fileparts(filepath), fileparts(unzipped_files{1}));

	% Remove the NUFFT root if it's in the path. Messing around with MEX files
	% that are potentially loaded can cause trouble.
	paths = strsplit(path, pathsep);
	if any(strcmp(nufft_root, paths))
		rmpath(nufft_root);
	end

	if force_compile
		% NOTE: It's important to run this here. If we delete MEX files that
		% are in the load path, Octave will crash.
		fprintf('Removing shipped MEX files.\n');
		delete(fullfile(nufft_root, '*.mex*'));
	end

	cd(nufft_root);
	try
		nufft1d1(1, 0, 0, 1, 1, 1);
		nufft2d1(1, 0, 0, 0, 1, 1, 1, 1);
		nufft3d1(1, 0, 0, 0, 0, 1, 1, 1, 1, 1);
		fprintf('MEX files are working.\n');
	catch
		fprintf('MEX files are not working. Compiling...\n')
		status = system('gfortran -fPIC -O2 -c *.f');
		if status == 1
			error('Compilation of Fortran files failed!');
		end

		mex nufft1d.c nufft1df90.o dirft1d.o dfftpack.o next235.o
		mex nufft2d.c nufft2df90.o dirft2d.o dfftpack.o next235.o
		mex nufft3d.c nufft3df90.o dirft3d.o dfftpack.o next235.o

		system('rm *.o');

		try
			nufft1d1(1, 0, 0, 1, 1, 1);
			nufft2d1(1, 0, 0, 0, 1, 1, 1, 1);
			nufft3d1(1, 0, 0, 0, 0, 1, 1, 1, 1, 1);
			fprintf('MEX files are working.\n');
		catch
			error('MEX compilation failed.');
		end
	end
	cd(current_dir);

	addpath(nufft_root);

	fprintf('CIMS NUFFT successfully installed!\n');
end
