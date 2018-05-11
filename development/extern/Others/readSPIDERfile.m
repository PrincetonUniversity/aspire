function data=readSPIDERfile(filename)
% readSPIDERfile : read a file in SPIDER format. It returns: 
%   a 2D matrix if input is a SPIDER image,
%   a 3D matrix if input is a SPIDER volume or a stack of images.
% 
% Returned matrix is type double. Does not handle complex formats.
%
% requires readSPIDERheader.m
%
% version 1.0 (Feb 2009) B. Baxter
% Copyright (C) 2009 Health Research Inc.
% Tested in Matlab 7.7.0 (R2008b, R2008a)

data = -1;

h = readSPIDERheader(filename);
if h == -1
  disp(sprintf('%s does not seem to be a valid SPIDER file\n', filename));
  return
end

nslice    = h(1);
nrow      = h(2);
nsam      = h(12);
nhdrbytes = h(22);  % aka labbyt
isstack   = h(24);
nstackitems = h(26);
if h(29) == 1       %  h(29):  1 = big-endian, 0 = little-endian
    endian = 'ieee-be';
else
    endian = 'ieee-le';
end

displayvals = 0;
if displayvals == 1
    disp(sprintf('nslice: %d', nslice));
    disp(sprintf('nrow: %d', nrow));
    disp(sprintf('nsam: %d', nsam));
    disp(sprintf('nhdrbytes: %d', nhdrbytes));
    disp(sprintf('isstack: %d', isstack));
end

[fp, errmsg] = fopen(filename, 'r', endian);

if fp == -1
  disp(errmsg)
  return
end

% read the header
nwords = nhdrbytes/4;
[hdr, count] = fread(fp, nwords, 'float32');

% SPIDER image file
if (nslice == 1) && (isstack == 0) 
    tdata = fread(fp, [nsam,nrow], 'float32');
    data = tdata';   % transpose the data
    fclose(fp);
    
% SPIDER volume    
elseif (nslice > 1) && (isstack == 0)
    % switch nsam and nrow since we're going to transpose the matrix
    data = zeros([nrow nsam nslice], 'single');
    for v = 1:nslice
       slice = fread(fp, [nsam,nrow], 'float32');
       data(:,:,v) = slice';
    end
    fclose(fp);
    
% stack of images    
elseif (nslice == 1) && (isstack > 0)
    % switch nsam and nrow since we're going to transpose the matrix
    data = zeros([nrow nsam nstackitems], 'single');
    for v = 1:nstackitems
       hdr = fread(fp, nwords, 'float32');
       slice = fread(fp, [nsam,nrow], 'float32');
       data(:,:,v) = slice';
    end
    fclose(fp);
end
