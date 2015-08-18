function cryo_workflow_preprocess_validate(workflow_fname)
% Validate that the workflow file has all the required parameters to
% exectue preprocessing.
%
% Yoel Shkolnisky August 2015.

% Read workflow file
tree=xmltree(workflow_fname);
workflow=convert(tree);

% Validate struct
assertfield(workflow,'info','working_dir');
assertfield(workflow,'info','logfile');
assertfield(workflow,'info','rawdata');

assertfield(workflow,'preprocess','phaseflip');
assertfield(workflow,'preprocess','ctfdata');
assertfield(workflow,'preprocess','nprojs');
assertfield(workflow,'preprocess','do_crop');
assertfield(workflow,'preprocess','croppeddim');
assertfield(workflow,'preprocess','do_downsample');
assertfield(workflow,'preprocess','downsampleddim');
assertfield(workflow,'preprocess','do_normalize');
assertfield(workflow,'preprocess','do_prewhiten');
assertfield(workflow,'preprocess','split');
assertfield(workflow,'preprocess','numgroups');