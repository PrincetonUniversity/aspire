function writeSPIDERfile(filename, data, varargin)
% writeSPIDERfile writes data to a file in SPIDER format.
%
% 2D data are written as a SPIDER image.
% 3D data may be written as either a SPIDER volume (default) or
% as a stack of 2D images (requires optional 3rd argument='stack').
%
% SPIDER files do not have a specific file extension.
% 
% Examples
% writeSPIDERfile('img001.dat', I)
% writeSPIDERfile('abc001.stk', V, 'stack')
%
% version 1.0 (Feb 2009) B. Baxter
% Copyright (C) 2009 Health Research Inc.
% Tested in Matlab 7.7.0 (R2008b, R2008a)

if (nargin > 2) && strcmp(varargin,'stack')
    stackarg = 0;      % make a stack 
    write_stack = 1;
else
    stackarg = -1;     % not a stack
    write_stack = 0;
end

datasize = size(data);

header = makeSPIDERheader(datasize, stackarg);
if header < 0
    disp('Unable to create header')
    return
end

[fp, errmsg] = fopen(filename, 'wb');

if fp == -1
  disp(errmsg)
  return
end

ndim = ndims(data);

if ndim == 2
    count = fwrite(fp, header, 'float32');
    count = fwrite(fp, data', 'float32');

elseif ndim == 3

   % write the volume header or overall stack header
   count = fwrite(fp, header, 'float32');
   
   nimages = datasize(3);
   nsam = datasize(1);
   nrow = datasize(2);
   imgsize = [nsam, nrow];
   
   for i = 1:nimages
       if write_stack
          hdr = makeSPIDERheader(imgsize, i);
          count = fwrite(fp, hdr, 'float32');
       end
       imgdata = data(:,:,i);
       count = fwrite(fp, imgdata', 'float32');
   end
   
else
    disp(sprintf('output data must be 2 or 3 dimensions\n'))
end
fclose(fp);

% -----------------------------------------------------

function header=makeSPIDERheader(datasize, stackarg)
% Create a header for a SPIDER file.
% 
% Stack files present a special case. They have an overall
% header, plus each image in the stack has its own header.
% In the overall header, istack > 0, maxim = total # of images. 
% In individual stack images, imgnum > 0
%
% stackarg: < 0 : not a stack
%           = 0 : make overall stack header
%           > 0 : imgnum in individual stack image
% Only image stacks are supported.

n = length(datasize);

nsam = datasize(2);   % switch due to row, column indexing
nrow = datasize(1);

if (n == 2) && (stackarg < 0)  % 2D image
    nslice = 1;
    iform = 1;
elseif (n == 3) && (stackarg < 0) % 3D volume
    nslice = datasize(3);
    iform = 3;
elseif stackarg == 0    % overall stack file header
    nslice = 1;
    iform = 1;
    istack = 2;
    maxim = datasize(3);
    imgnum = 0;
elseif stackarg > 0    % image within a stack file
    nslice = 1;
    iform = 1;
    istack = 0;
    maxim = 0;
    imgnum = stackarg;
end

lenbyt = nsam * 4;
labrec = uint32(1024 / lenbyt);
if mod(1024,lenbyt) ~= 0
    labrec = labrec + 1;
end
labbyt = labrec * lenbyt;
nvalues = uint32(labbyt/4);

if nvalues < 256
    header = -1;
    return
end

header = zeros([nvalues 1], 'single');

header(1)  = nslice;  % nslice (=1 for an image) 
header(2)  = nrow;    % number of rows per slice
header(5)  = iform;   % iform: 1 for 2D image, 3 for volume
header(12) = nsam;    % number of pixels per line
header(13) = labrec;  % number of records in file header
header(22) = labbyt;  % total number of bytes in header
header(23) = lenbyt;  % record length in bytes
if stackarg >= 0
   header(24) = istack;  % > 0 in overall stack file header
   header(26) = maxim;   % total number of images in a stack file
   header(27) = imgnum;  % current image in a stack file
end