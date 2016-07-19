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

%% Read preprocessing parameters
% Do phase flip?
message='Phaseflip projections? ';
do_phaseflip=multichoice_question(message,{'Y','N'},[ 1, 0],'Y');
if do_phaseflip==1
    ctfdata=fmtinput('Enter full path of CTF file: ','','%s');
end

% Do downsample?
[~,mrc_header]=ReadMRC(workflow.info.rawdata,1, 1);
% Input file has  mrc_header.nz images, each of size mrc_header.nx x
% mrc_header.ny
fprintf('%s :%dx%d (%d projections)\n',...
    workflow.info.rawdata,mrc_header.nx,mrc_header.ny,mrc_header.nz);
message='Number of projections to read? ';
nprojs=fmtinput(message,mrc_header.nz,'%d');

croppeddim=mrc_header.nx; % Default size of cropped projections is the current size.
message='Crop? ';
do_crop=multichoice_question(message,{'Y','N'},[ 1, 0],'Y');
if do_crop==1
    message='Crop to size? ';
    croppeddim=fmtinput(message,mrc_header.nx,'%d');
end


downsampleddim=mrc_header.nx; % Default size of downsampled projections is the current size.
message='Downsample? ';
do_downsample=multichoice_question(message,{'Y','N'},[ 1, 0],'Y');
if do_downsample==1
    message='Downsample to size? ';
    downsampleddim=fmtinput(message,mrc_header.nx,'%d');
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
workflow.preprocess.ctfdata='';
if do_phaseflip
    workflow.preprocess.ctfdata=ctfdata;
end
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
