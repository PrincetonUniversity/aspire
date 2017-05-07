function star2mrc(starname,mrcname,basedir,ignorepath)
%
% STAR2MRC     Create image stack from STAR file
%
% star2mrc(starname,mrcname,basedir)
%   Create an MRC file name mrcname from images specified in the STAR file
%   starname. The relative paths of the images in the STAR file relative to
%   basedir. Set ignorepath to nonzero to ignore the image path specified
%   in the STAR file. In such a case only the image name is used, relative
%   to basedir. This option is useful if the stack files have been moved to
%   a directory different than the one specified in the STAR file.
%
% Example:
%   star2mrc('shiny_2sets.star','rawdata.mrc','~/datasets/100028/');
%
% Yoel Shkolnisky, December 2016.

if ~exist('mrcname','var')
    error('Name of MRC output stack not given');
end

if ~exist('basedir','var')
    basedir='.';
end

if ~exist('ignorepath','var')
    ignorepath=0;
end

% Load the information from the star file.
log_message('Loading %s',starname);
datablocks=readSTAR(starname);

N=numel(datablocks.data); % Number of images to read
stackwriter=imagestackWriter(mrcname,N,1,1000);
prev_micrographname=[]; % Name of previously loaded micrograph

% Read images appearing in the star file
log_message('Copying images to MRC %s',mrcname);
printProgressBarHeader;
for i=1:N
    progressTicFor(i,N);
    imkey=datablocks.data{i}.rlnImageName; % String encoding micrograph 
                          % name and image index withing the micrograph
    keyparts=strsplit(imkey,'@');
    idx=str2double(keyparts{1});
    micrographname=keyparts{2};
   
    if ignorepath
        [~,micrographname,ext]=fileparts(micrographname);
        micrographname=[micrographname ext];
    end
    
    if ~strcmp(micrographname,prev_micrographname)
        prev_micrographname=micrographname;
        fname=fullfile(basedir,micrographname);
        micrographreader=imagestackReader(fname,1000);
    end
                    
    im=micrographreader.getImage(idx);
    if mod(size(im,1),2)==0
       im=im(1:end-1,1:end-1);
    end	
    stackwriter.append(im);
end
stackwriter.close;
