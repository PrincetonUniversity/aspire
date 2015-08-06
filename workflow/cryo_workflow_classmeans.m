function cryo_workflow_classmeans(workflow_fname)
% CRYO_WORKFLOW_CLASSMEANS  Interactive generation of class means
%
% cryo_workflow_classify(workflow_fname)
%   Interactively collect all parameters required to generate class means
%   from precomputed classification data, and generate the class means.
%   workflow_name is the name of the file where the entered parameters will
%   be saved.
%
% See also cryo_workflow_start
%
% Yoel Shkolnisky, August 2015.

% Read workflow file
tree=xmltree(workflow_fname);
workflow=convert(tree);

%% Read Classification parameters

maxavg=str2double(workflow.classification.k_VDM_out);
defnnavg=50;
message=sprintf('Number of projections to average into each class means (up to %d)? ',maxavg);
nnavg=fmtinput(message,defnnavg,'%d');

if nnavg>maxavg
    log_message('Cannot average %s projections. Max averaging allowed is %d.',nnavg,maxavg);
    log_message('Run cryo_workflow_classify again and increase k_VDM_out to at least %d',nnavg);
    log_message('Aborting...');
    return;
end

%% Update workflow struct

workflow.classmeans.nnavg = nnavg; % output number of nearest neighbors

tree=struct2xml(workflow);
save(tree,workflow_fname); 

%% Generate class means
cryo_workflow_classmeans_execute(workflow_fname);
