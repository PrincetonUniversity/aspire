function cryo_workflow_classmeans_validate(workflow_fname)
% Validate that the section 'preprocess' in the given XML workflow has all
% the required fields to exectue preprocessing.
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
