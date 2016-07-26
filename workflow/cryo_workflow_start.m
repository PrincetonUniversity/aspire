function cryo_workflow_start
%
% CRYO_WORKFLOW_START   Start ASPIRE processing workflow
%
% cryo_workflow_start
%   Gather all information required to start processing workflow of raw
%   images.
% This is the first command to be called in any processing workflw.
% 
% Yoel Shkolnisky, August 2015.

workflow_name='';
while isempty(workflow_name)
    workflow_name=fmtinput('Enter workflow name: ','','%s');
    if ~isvarname(workflow_name)
        fprintf('Workflow name must be a valid filename.\n');
        workflow_name='';
    end
end
workflow_desc=fmtinput('Enter workflow description: ','','%s');

workflow_mrc='';
while isempty(workflow_mrc)
    workflow_mrc =fmtinput('Enter full path input MRC file: ','','%s');
    if exist(workflow_mrc,'file')~=2
        fprintf('MRC file does not exist.\n');
        workflow_mrc='';
    end
end

workflow_dir =fmtinput('Enter full path of output directory: ','','%s');
if ~exist(workflow_dir,'dir') % Do we need to create directory?
    message='Output directory does not esxist. Create?';
    do_create=multichoice_question(message,{'Y','N'},[ 1, 0],'Y');
    if do_create==1
        mkdir(workflow_dir);
    end
else % Check if directory empty
    listing=dir(workflow_dir);
    if numel(listing)>2
        message='Output directory not empty. Continue?';
        do_continue=multichoice_question(message,{'Y','N'},[ 1, 0],'N');
        if do_continue==0
            fprintf('Aborting...\n');
            return;
        end
    end
end

workflow_log =fmtinput('Enter log file name: ','log.txt','%s');

fprintf('Preparing data. Please wait...\n');


% Next lines have been commented until we decided we want to keep the hash
% of the input data set. This requires finding a freely available hash
% function which is fast enough for large inputs. 
% Y.S. July 2016.

% % % % Compute hash of data file
% % % opt.Input='file';
% % % opt.Format='hex';
% % % opt.Method='MD5';
% % % hash=DataHash(workflow_mrc,opt);

% Create strucut
workflow.info.name=workflow_name;
workflow.info.description=workflow_desc;
workflow.info.created=datestr(now);
workflow.info.rawdata=workflow_mrc;
% % % workflow.info.rawdatahash=hash;
workflow.info.working_dir=workflow_dir;
workflow.info.logfile=workflow_log;

% Save as xml
tree=struct2xml(workflow);
xmlname=sprintf('%s.xml',workflow_name);
fname=fullfile(workflow_dir,xmlname);
save(tree,fname); 

fprintf('Workflow file: %s\n',fname);
fprintf('Use this file name when calling subsequent funtions.\n');
fprintf('Call next cryo_workflow_preprocess(''%s'')\n',fname);
