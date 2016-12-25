function star2mrc(starname,mrcname,basedir)
%
% STAR2MRC     Create image stack from STAR file
%
% star2mrc(starname,mrcname,basedir)
%   Create an MRC file name mrcname from images specified in the STAR file
%   starname. The relative paths of the images in the STAR file relative to
%   basedir.
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

% Load the information from the star file.
log_message('Loading %s',starname);
datablocks=readSTAR(starname);

N=numel(datablocks.data); % Number of images to read
stackwriter=imagestackWriter(mrcname,1000,N);
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
    
    if ~strcmp(micrographname,prev_micrographname)
        prev_micrographname=micrographname;
        fname=fullfile(basedir,micrographname);
        micrographreader=imagestackReader(fname,1000);
    end
                    
    im=micrographreader.getImage(idx);
    stackwriter.append(im);
end
stackwriter.close;