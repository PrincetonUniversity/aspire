% INSTALL_FINUFFT Install the Flatiron Institute NUFFT library
%
% Usage
%    install_finufft();
%    install_finufft(url, location, fftw_location);
%
% Input
%    url: The url from which the package should be downloaded. By default,
%       this is
%
%          https://github.com/ahbarnett/finufft/archive/
%             a960a0f17e196864e6f6bede8bdc21bb8ff5e653.zip
%
%       Which is corresponds to version 0.9 of the library released on
%       2017-06-17.
%    location: The location in which the package should be installed. By
%       default, this is the subfolder 'extern' of the ASPIRE root folder.
%    fftw_location: The location of the FFTW3 package needed to compile FINUFFT.
%       By default, the location is 'extern/fftw3'.
%
% Description
%    This function downloads and installs the FINUFFT library from the Flatiron
%    Institute, which implements a non-uniform fast Fourier transform.

function install_finufft(url, location, fftw_location)
	if nargin < 1 || isempty(url)
		url = ['https://github.com/ahbarnett/finufft/archive/' ...
			   'a960a0f17e196864e6f6bede8bdc21bb8ff5e653.zip'];
	end

	if nargin < 2 || isempty(location)
		location = fullfile(aspire_root(), 'extern');
	end

	if nargin < 3 || isempty(fftw_location)
		fftw_location = fullfile(aspire_root(), 'extern', 'fftw3');
	end

	if ~exist(fftw_location, 'dir')
		fprintf('FFTW3 not installed. Installing...\n');
		install_fftw();
	end

	fprintf('Installing the Flatiron Institute NUFFT package.\n');
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

	% Can't use output from unzip since it seems to contain non-existent files?
	% Probably an Octave bug.
	[~, commit_id, ~] = fileparts(filepath);
	finufft_dir = ['finufft-' commit_id];

	old_finufft_root = fullfile(fileparts(filepath), finufft_dir);
	new_finufft_root = fullfile(fileparts(filepath), 'finufft');

	if exist(new_finufft_root, 'dir')
		error(['Unable to move folder into ''extern/finufft''. Perhaps ' ...
			'FINUFFT is' char(10) 'already installed? If you wish to ' ...
			'reinstall, please remove this folder.']);
	end

	movefile(old_finufft_root, new_finufft_root);
	finufft_root = new_finufft_root;

	cd(finufft_root);

	cmd = 'make';

	if ~isoctave()
		cmd = [cmd ' matlab'];
	else
		cmd = [cmd ' octave'];
	end

	old_cpath = getenv('CPATH');
	new_cpath = [fullfile(fftw_location, 'include') pathsep ...
		old_cpath];
	setenv('CPATH', new_cpath);

	old_library_path = getenv('LIBRARY_PATH');
	new_library_path = [fullfile(fftw_location, 'lib') pathsep ...
		old_library_path];
	setenv('LIBRARY_PATH', new_library_path);

	status = system(cmd);

	setenv('CPATH', old_cpath);
	setenv('LIBRARY_PATH', old_library_path);

	if status ~= 0
		error('''make'' failed');
	end

	cd(current_dir);

	addpath(fullfile(location, 'finufft', 'matlab'));
	try
		finufft1d1(0, 0, 1, 1, 1);
	catch
		error('Flatiron Institute NUFFT installation failed\n');
		rmpath(fullfile(location, 'finufft', 'matlab'));
	end

	fprintf('Flatiron Institute NUFFT successfully installed!\n');
end
