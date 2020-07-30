function star2mrc(starname,mrcname,basedir,ignorepath,Nprojs,verbose)
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
%   verbose=0   silent
%   verbose=1   show progress bar
%   verbose=2   print progress every 1000 images
%   verbose=3   print message for each processed image
%
% Example:
%   star2mrc('shiny_2sets.star','rawdata.mrcs','~/datasets/100028/');
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

if ~exist('Nprojs','var')
    Nprojs=-1;
end

if ~exist('verbose','var')
    verbose=0;
end

% Load the information from the star file.
log_message('Loading %s',starname);
datablocks=readSTAR(starname);

% NOTE: If STAR files of RELION 3.1 is used, then the structure of the
% STAR file is assumed to contained one optics group (location 1 in the
% stardata array) and one particles group (location 2 in the stardata
% array).
if numel(datablocks)==1 % RELION version < 3.1
    Nrecords=numel(datablocks.data);
else
    Nrecords=numel(datablocks(2).data);
end


if Nprojs~=-1
    N=min(Nprojs,Nrecords); % Number of images to read
else
    N=Nrecords;
end

if verbose>=1
    log_message('Reading %d projections',N);
end

stackwriter=imagestackWriter(mrcname,N,1,1000);
prev_stackname=[]; % Name of previously loaded micrograph

% Read images appearing in the star file
if verbose>0
    log_message('Copying images to MRC %s',mrcname);
end

if verbose==1
    printProgressBarHeader;
end

for i=1:N
    if verbose==1
        progressTicFor(i,N);
    end
        
    if numel(datablocks)==1 % RELION version < 3.1        
        imkey=datablocks.data{i}.rlnImageName; % String encoding micrograph
        % name and image index withing the micrograph
    else
        imkey=datablocks(2).data{i}.rlnImageName;
    end
    
    keyparts=strsplit(imkey,'@');
    idx=str2double(keyparts{1});
    stackname=keyparts{2};

    if verbose>=3
        log_message('Processing prjection %d/%d: reading image %d from stack %s',...
            i,N,idx,stackname);
    end
    
    
    if ignorepath
        [~,stackname,ext]=fileparts(stackname);
        tmp=[stackname ext];
        stackname=tmp;
    end
    
    if ~strcmp(stackname,prev_stackname)
        prev_stackname=stackname;
        fname=fullfile(basedir,stackname);
        micrographreader=imagestackReader(fname,1000);
        
        if verbose>=2
            log_message('Processing prjection %d/%d: Initializing MRC stack %s',...
                i,N,fname);
        end
        
    end
                    
    im=micrographreader.getImage(idx);
    if mod(size(im,1),2)==0
       im=im(1:end-1,1:end-1);
    end	
    stackwriter.append(im);
    
    if verbose==2
        if mod(i,1000)==0
            log_message('Finsihed processing %d/%d projections',i,N);
        end
    end
    
end
stackwriter.close;
