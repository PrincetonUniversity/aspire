function cryo_workflow_classmeans_validate(workflow_fname)
% Validate that the workflow file has all the required parameters to
% generate class means.
%
% Yoel Shkolnisky August 2015.

% Read workflow file
tree=xmltree(workflow_fname);
workflow=convert(tree);

% Validate struct

assertfield(workflow,'info','working_dir');
assertfield(workflow,'info','logfile');
assertfield(workflow,'preprocess','numgroups');
assertfield(workflow,'classmeans','nnavg');
assertfield(workflow,'classmeans','num_averages');
assertfield(workflow,'classmeans','use_EM');

