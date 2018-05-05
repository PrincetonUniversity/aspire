% PICK_NUFFT_LIBRARY Specifies which NUFFT library to use for a given size
%
% Usage
%     lib_code = pick_nufft_library(sz);
%
% Input
%     sz: The size of the Fourier transform box.
%
% Output
%     lib_code: The library to use.
%
% Description
%    This is an internal function that should not be used directly. Please
%    see get_nufft_library and set_nufft_library.
%
% See also
%    have_nufft_library

function lib_code = pick_nufft_library(sz)
	lib_codes = nufft_library_codes();

	lib_codes = [lib_codes 1];

	if numel(sz) == 2 && sz(2) == 1
		sz = sz(1);
	end

	if ~all(mod(sz, 2)==0) || any(sz<10)
		lib_codes(lib_codes==3) = [];
	end

	while ~have_nufft_library(lib_codes(1))
		lib_codes(1) = [];
	end

	lib_code = lib_codes(1);

	if lib_code == 1
		warning('aspire:using_nudft', ...
			['Using naive DFT implementation for non-uniform discrete ' ...
			'Fourier transform. Please install an NUFFT library by running ' ...
			'''install_cims_nufft'' or ''install_chemnitz_nfft''.']);
	end
end
