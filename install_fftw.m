% INSTALL_FFTW Install FFTW library
%
% Usage
%    install_fftw();
%    install_fftw(url, location);
%
% Input
%    url: The url from which the package should be downloaded. By default,
%       this is
%
%          ftp://ftp.fftw.org/pub/fftw/fftw-3.3.4.tar.gz
%
%    location: The location in which the package should be installed. By
%       default, this is the subfolder 'extern' of the ASPIRE root folder.
%
% Description
%    This function downloads and compiles the FFTW library, which implements
%    fast Fourier transforms.

function install_fftw(url, location)
	if nargin < 1 || isempty(url)
		url = 'ftp://ftp.fftw.org/pub/fftw/fftw-3.3.4.tar.gz';
	end

	if nargin < 2 || isempty(location)
		location = fullfile(aspire_root(), 'extern');
	end

	fprintf('Installing the FFTW package.\n');
	fprintf('URL: %s\n', url);
	fprintf('Location: %s\n', location);
	fprintf('\n');

	ind = find(url=='/', 1, 'last');
	filename = url(ind+1:end);

	filepath = fullfile(location, filename);

	if exist(filepath, 'file')
		fprintf('Package exists on disk. Skipping download.\n');
	else
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

	fftw_dir = fileparts(unzipped_files{1});
	ind = find(fftw_dir=='/', 1);
	if ~isempty(ind)
		fftw_dir = fftw_dir(1:ind-1);
	end

	fftw_root = fullfile(fileparts(filepath), fftw_dir);

	cd(fftw_root);

	status = system(['./configure ' ...
	                 '--prefix=' fullfile(location, 'fftw3') ' ' ...
	                 '--enable-openmp --enable-shared --enable-threads']);
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

	fprintf('FFTW successfully installed!\n');
end
