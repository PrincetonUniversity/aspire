% HAVE_NUFFT_LIBRARY Determines whether an NUFFT library is installed
%
% Usage
%    have_library = have_nufft_library(lib_code);
%
% Input
%    lib_code: A library code. See documentation of nufft_library_code for
%       more information.
%
% Output
%    have_library: True if the specified library is installed.
%
% Description
%    This is an internal function that should not be used directly. Please
%    see get_nufft_library and set_nufft_library.
%
% See also
%    nufft_library_code

function have_library = have_nufft_library(lib_code)
	if lib_code == 1
		have_library = true;
	elseif lib_code == 2
		have_library = true;
		try
			nufft1d1(1, 0, 0, 1, 1, 1);
		catch
			have_library = false;
		end
	elseif lib_code == 3
		have_library = true;
		try
			p = nfft_init_1d(2, 1);
			nfft_finalize(p);
		catch
			have_library = false;
		end
	elseif lib_code == 4
		have_library = true;
		try
			finufft1d1(0, 0, 1, 1, 1);
		catch
			have_library = false;
		end
	else
		error('invalid library code');
	end
end
