function [m totalImages]=BufferedReadImage(fileName,index,chunk)
% function [m totalImages]=BufferedReadImage(filename,index,chunk)
% Read a single image from an Imagic or MRC file, while keeping a cache to
% allow efficient loading of subsequent images too.  The file type is
% deduced from the extension (.hed, .img, .mrc).
% - The static buffer is deallocated if you call this function with no index
% argument, or if index <1.
% - The optional argument chunk gives the number of images to hold in buffer
% (default is 1000).
% As the OS does a good job of buffering reads, the local buffering is
% not strictly necessary, but it avoids overhead in our ReadImagic and 
% ReadMRC functions.

persistent imgs nTotal startImage nChunk name;

defaultChunk=1000;
totalImages=0;
m=[];

% No index or negative index means the buffer is cleared.
if nargin<2 || index<1
        imgs=[];
        return
end;

% If the filename changes, clear the buffer and store the name.
if ~strcmp(fileName,name)
    name=fileName;
    imgs=[];
end;

% If the chunk value is given, store it.
if nargin>2
    nChunk=chunk;
end;

nim=size(imgs,3); % Get the number of stored images (=1 if nim is empty)
if numel(imgs)<1 || index>startImage+nim-1 || index<startImage  % need data
    if numel(nChunk)<1  % hasn't been initialized
        nChunk=defaultChunk;
    end;
    [pa nm ext]=fileparts(fileName);

    % Read data starting at index
    switch lower(ext)
        case '.mrc'
            [imgs s]=ReadMRC(fileName,index,nChunk);
            startImage=index;
            nTotal=s.nz;
        case {'.hed','.img'}
            [imgs info]=ReadImagic(fileName,index,nChunk);
            startImage=index;
            nTotal=info.nim;
        otherwise
            error(['Invalid file extension in ' fileName]);
    end;
end;
totalImages=nTotal;  % copy the returned value
copyIndex=index-startImage+1;
if copyIndex <= size(imgs,3)
    m=imgs(:,:,index-startImage+1); % get an image from the buffer
end;