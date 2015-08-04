function cryo_workflow_classmeans(workflow_fname)

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
