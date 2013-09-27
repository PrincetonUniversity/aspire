function [m pixA ok]=ReadEMFile(name)
% function [m pixA ok]=ReadEMFile(name)
% General image file reader.
% Based on the file extension, read the image and return the pixel size
% (in angstroms).  The returned value of pixA is zero if a pixel size is
% not available.
% Note that the original data type (uint8, uint16, single) is returned!  No
% conversions are done.  However, rotations are made to standard graphics
% files to correspond to the Cartesian coordinate convention.
% This function understands the following extensions:
% .tif  (generic, or TVIPS files) --rotated 270 degrees
% .dm3
% .mrc
% .hed, .img
% .jpg   --rotated 270 degrees
% If the file has the wrong extension, m is returned as zeros and ok=0.
% fs Jan 2011

[pa nm ex]=fileparts(name);
pixA=0;  % default value
m=zeros(100,100); % default value
ok=1;
switch lower(ex)
    case '.mrc'
        [m s]=ReadMRC(name);
        pixA=s.rez/s.nx;
    case '.dm3'
        [m pixnm]=ReadDM3(name);
        pixA=10*pixnm;
    case '.tif'
        m=rot90(imread(name),3);
        s=ReadTiffMetadata(name);  % Get metadata from a TVIPS file.
        pixA=s.pixA;               % will be zero if not a TVIPs file
%         Although TIFF doesn't support it, TVIPS 2-byte files are encoded
%         as signed integers, so we coerce the type.
        if isa(m,'uint16')&&(pixA>0)
            sz=size(m);
            m=typecast(m(:),'int16');
            m=reshape(m,sz);
        end;
    case  {'.hed','.img'}
        m=ReadImagic(name);
    case '.jpg'
        m=rot90(imread(name),3);  % allow jpegs to be read too.
    otherwise
        warning(['Unknown file type: ' ex]);
        ok=0;
end;
