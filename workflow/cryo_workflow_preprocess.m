function cryo_workflow_preprocess(workflow_fname)
% CRYO_WORKFLOW_PREPROCESS  Interactive data set preprocessing
%
% cryo_workflow_preprocess(workflow_fname)
%   Interactively collect all parameters required to preprocess the
%   projections data set and execute the preprocessing.
%   workflow_name is the name of the file where the entered parameters will
%   be saved.
%
% See also cryo_workflow_start
%
% Yoel Shkolnisky, August 2015.

% Read workflow file
tree=xmltree(workflow_fname);
workflow=convert(tree);

% Check that the input has not changed since the last time is was
% preprocessed. % Compute hash of data file and compare to the stored hash.
%fprintf('Validating input...\n');
%opt.Input='file';
%opt.Format='hex';
%opt.Method='MD5';
% hash=DataHash(workflow.info.rawdata,opt);
% if ~strcmpi(hash,workflow.info.rawdatahash)
%     message='Raw projection''s file has changed. Proceed? ';
%     do_proceed=multichoice_question(message,{'Y','N'},[ 1, 0],'Y');
%     if ~do_proceed
%        log_message('Aborting...\n');
%        return;
%     end
% end

%[data_dir,~,data_ext]=fileparts(workflow.info.rawdata);
[~,~,data_ext]=fileparts(workflow.info.rawdata);
data_dir=[]; % Look for stack files exatcly as appears in the STAR file.

%% Read preprocessing parameters

do_phaseflip=0;
pixA=-1;
if ~strcmpi(data_ext,'.star')
    fprintf('Input data is not a STAR file. Skipping phase flipping\n');
else        
    % Do phase flip?
    message='Phaseflip projections? ';
    do_phaseflip=multichoice_question(message,{'Y','N'},[ 1, 0],'Y');
    % Read pixel size
    if do_phaseflip
        pixA=fmtinput('Enter pixel size in Angstrom (-1 to read from STAR file): ',-1,'%f');
    end
end

if strcmpi(data_ext,'.star')
    fprintf('Reading STAR file %s...\n',workflow.info.rawdata);
    stardata=readSTAR(workflow.info.rawdata,1);
    
    % Read the first image from the star data to read the dimensions of the
    % images.
    % NOTE: If STAR files of RELION 3.1 is used, then the structure of the
    % STAR file is assumed to contained one optics group (location 1 in the
    % stardata array) and one particles group (location 2 in the stardata
    % array). 
    if numel(stardata)==1 % RELION version < 3.1
        nz=numel(stardata.data);
        imageID=stardata.data.rlnImageName;
    else % RELION 3.1
        nz=numel(stardata(2).data);
        imageID=stardata(2).data{1}.rlnImageName;
    end
    imparts=strsplit(imageID,'@');
    %imageidx=str2double(imparts{1});
    stackname=imparts{2};
    
    MRCname=fullfile(data_dir,stackname);
    projs1=ReadMRC(MRCname,1,1);
    nx=size(projs1,1);
    ny=size(projs1,2);
else
    fprintf('Reading MRC file %s...\n',workflow.info.rawdata);
    [~,mrc_header]=ReadMRC(workflow.info.rawdata,1, 1);
    % Input file has  mrc_header.nz images, each of size mrc_header.nx x
    % mrc_header.ny
    nx=mrc_header.nx;
    ny=mrc_header.ny;
    nz=mrc_header.nz;
end

if nx~=ny
    error('Input projections must be square.')
end


% Do downsample?
fprintf('%s :%dx%d (%d projections)\n',...
    workflow.info.rawdata,nx,ny,nz);
message='Number of projections to read? ';
nprojs=fmtinput(message,nz,'%d');

croppeddim=nx; % Default size of cropped projections is the current size.
message='Crop? ';
do_crop=multichoice_question(message,{'Y','N'},[ 1, 0],'Y');
if do_crop==1
    message='Crop to size? ';
    croppeddim=fmtinput(message,nx,'%d');
end


downsampleddim=nx; % Default size of downsampled projections is the current size.
message='Downsample? ';
do_downsample=multichoice_question(message,{'Y','N'},[ 1, 0],'Y');
if do_downsample==1
    message='Downsample to size? ';
    downsampleddim=fmtinput(message,nx,'%d');
end

% Do normalize background.
message='Normalize background of images to variance 1? ';
do_normalize=multichoice_question(message,{'Y','N'},[ 1, 0],'Y');

% Do prewhiten
message='Prewhiten? ';
do_prewhiten=multichoice_question(message,{'Y','N'},[ 1, 0],'Y');

% Do split into groups?
do_shuffle=0;
message='Split data into groups? ';
do_split=multichoice_question(message,{'Y','N'},[ 1, 0],'Y');
if do_split==1
    numgroups=fmtinput('Number of groups ',2,'%d');
    
    % Split data at random?
    message='Shuffle data before splitting? ';
    do_shuffle=multichoice_question(message,{'Y','N'},[ 1, 0],'N');
else
    numgroups=1;
end

%% Update workflow struct
workflow.preprocess.phaseflip=do_phaseflip;
workflow.preprocess.pixA=pixA;
workflow.preprocess.nprojs=nprojs;
workflow.preprocess.do_crop=do_crop;
workflow.preprocess.croppeddim=croppeddim;
workflow.preprocess.do_downsample=do_downsample;
workflow.preprocess.downsampleddim=downsampleddim;
workflow.preprocess.do_normalize=do_normalize;
workflow.preprocess.do_prewhiten=do_prewhiten;
workflow.preprocess.split=do_split;
workflow.preprocess.numgroups=numgroups;
workflow.preprocess.do_shuffle=do_shuffle;

tree=struct2xml(workflow);
save(tree,workflow_fname); 

%% Execute preprocessing
cryo_workflow_preprocess_execute(workflow_fname);
