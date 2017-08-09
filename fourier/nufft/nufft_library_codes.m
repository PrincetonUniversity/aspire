% NUFFT_LIBRARY_CODES Get/set NUFFT library codes
%
% Usage
%    lib_codes = nufft_library_codes();
%    old_lib_codes = nufft_library_codes(new_lib_codes);
%
% Input
%    new_lib_codes: The new library code order.
%
% Output
%    old_lib_codes: The old library code order.
%
% Description
%    This is an internal function that should not be used directly. Please
%    see get_nufft_library and set_nufft_library.
%
%    The library code order specifies in which order to try the different
%    NUFFT libraries. The code values are:
%       1: The naive DFT implementation, 'dft',
%       2: The CIMS NUFFT library by Greengard and Lee, 'cims', and
%       3: The TU Chemnitz NFFT library by Potts, 'chemnitz'.
%       4: The Flatiron Institure NUFFT library by Barnett and Magland,
%          'finufft'.
%    For example, the default library code order, [3 2 1] specifies to first
%    try the 'chemnitz' library, then 'cims', and if none of these work,
%    'nudft'.
%
% See also
%    get_nufft_libraries, set_nufft_libraries

function old_lib_codes = nufft_library_codes(new_lib_codes)
	persistent lib_codes;

	if isempty(lib_codes)
		lib_codes = [3 2 1];
	end

	old_lib_codes = lib_codes;

	if nargin == 1
		lib_codes = new_lib_codes;
	end
end
