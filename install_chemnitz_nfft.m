% INSTALL_CHEMNITZ_NFFT Install NFFT library from TU Chemnitz
%
% Usage
%    install_chemnitz_nfft();
%    install_chemnitz_nfft(url, location, fftw_location);
%
% Input
%    url: The url from which the package should be downloaded. By default,
%       this is
%
%          https://www-user.tu-chemnitz.de/~potts/nfft/download/
%             nfft-3.3.1.tar.gz
%
%    location: The location in which the package should be installed. By
%       default, this is the subfolder 'extern' of the ASPIRE root folder.
%    fftw_location: The location of the FFTW3 package needed to compile NFFT.
%       By default, the location is 'extern/fftw3'.
%
% Description
%    This function downloads and compiles the NFFT library from TU Chemnitz,
%    which implements a non-uniform fast Fourier transform.

function install_chemnitz_nfft(url, location, fftw_location)
	if nargin < 1 || isempty(url)
		url = ['https://www-user.tu-chemnitz.de/~potts/nfft/download/' ...
		       'nfft-3.3.2.tar.gz'];
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

	fprintf('Installing the TU Chemnitz NFFT package.\n');
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
	if numel(unzipped_files) == 1 && ...
	   strcmp(unzipped_files{1}(end-3:end), '.tar')
		untarred_files = untar(unzipped_files{1});
		delete(unzipped_files{1});
		unzipped_files = untarred_files;
	end

	nfft_dir = fileparts(unzipped_files{1});
	ind = find(nfft_dir=='/', 1);
	if ~isempty(ind)
		nfft_dir = nfft_dir(1:ind-1);
	end

	nfft_root = fullfile(fileparts(filepath), nfft_dir);

	cd(nfft_root);

	cmd = ['./configure ' ...
	         '--prefix=' fullfile(location, 'nfft') ' ' ...
	         '--enable-openmp --disable-applications ' ...
	         '--with-fftw3=' fftw_location ' '];

	if ~isoctave()
	    cmd = [cmd '--with-matlab=' matlabroot];
	else
	    cmd = [cmd '--with-octave=' matlabroot];
	end

	status = system(cmd);

	if status ~= 0
		error('''configure'' failed');
	end

	status = system('make');
	if status ~= 0
		error('''make'' failed\n');
	end

	status = system('make install');
	if status ~= 0
		error('''make install'' failed\n');
	end

	cd(current_dir);

	addpath(fullfile(location, 'nfft', 'lib'));
	addpath(fullfile(location, 'nfft', 'share', 'nfft', 'matlab', 'nfft'));
	try
		p = nfft_init_1d(2, 1);
		nfft_finalize(p);
	catch
		error('TU Chemnitz NFFT installation failed\n');
		rmpath(fullfile(location, 'nfft', 'lib'));
		rmpath(fullfile(location, 'nfft', 'share', 'nfft', 'matlab', 'nfft'));
	end

	fprintf('TU Chemnitz NFFT successfully installed!\n');
end
