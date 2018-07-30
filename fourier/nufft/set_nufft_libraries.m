% SET_NUFFT_LIBRARIES Set the order of NUFFT libraries to use
%
% Usage
%    set_nufft_libraries(new_libs);
%
% Input
%    new_libs: The new library order to use. This is a cell array containing
%       the following strings:
%          'dft': The naive DFT implementation,
%          'cims': The CIMS NUFFT implementation by Greengard & Lee, or
%          'chemnitz': The TU Chemnitz NFFT implementation by Potts.
%          'finufft': The Flatiron Institure NUFFT library by Barnett and
%             Magland.
%
% Description
%    The library order determines which libraries will be used when calling
%    the NUFFT wrappers. If one is not available, or is unable to process
%    a given array, the next one is tried.
%
% See also
%    get_nufft_libraries

function set_nufft_libraries(new_libs)
	lib_names = {'dft', 'cims', 'chemnitz', 'finufft'};

	if ~iscell(new_libs)
		new_libs = {new_libs};
	end

	for k = 1:numel(new_libs)
		new_lib_code(k) = find(strcmp(new_libs{k}, lib_names), 1);
	end

	nufft_library_codes(new_lib_code);
end
