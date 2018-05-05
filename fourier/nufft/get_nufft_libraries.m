% GET_NUFFT_LIBRARIES Get the order of NUFFT libraries to use
%
% Usage
%    libs = get_nufft_libraries();
%
% Output
%    libs: A cell array of the libraries that are to be used by the NUFFT
%       wrappers, in order of priority. A list of the compatible libraries
%       is given in the documentation for set_nufft_libraries.
%
% See also
%    set_nufft_libraries

function libs = get_nufft_libraries()
	lib_names = {'dft', 'cims', 'chemnitz', 'finufft'};

	libs = lib_names(nufft_library_codes());
end
