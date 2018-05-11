% GET_FILE_SIZE Obtain size of a given file
%
% Usage
%    bytes = get_file_size(filename);
%
% Input
%    filename: The name of a file whose size we would like to obtain.
%
% Output
%    bytes: The size of the file in bytes.

function bytes = get_file_size(filename)
    info = dir(filename);
    bytes = info.bytes;
end
