function cryo_workflow_classify_validate(workflow_fname)
% Validate that the workflow file has all the required parameters to
% exectue classification.
%
% Yoel Shkolnisky August 2015.

% Read workflow file
tree=xmltree(workflow_fname);
workflow=convert(tree);

% Validate struct
assertfield(workflow,'info','working_dir');
assertfield(workflow,'info','logfile');
assertfield(workflow,'preprocess','numgroups');
assertfield(workflow,'classification','n_nbor'); 
assertfield(workflow,'classification','isrann'); 
assertfield(workflow,'classification','k_VDM_in');
assertfield(workflow,'classification','VDM_flag');
assertfield(workflow,'classification','k_VDM_out');
