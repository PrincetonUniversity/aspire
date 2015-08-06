function cryo_workflow_abinitio(workflow_fname)
% CRYO_WORKFLOW_ABINITIO  Interactive generation of abinitio models
%
% cryo_workflow_classify(workflow_fname)
%   Interactively collect all parameters required to generate abinitio
%   models from precomputed class averages and execute the reconstruction.
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

defnmeans=1000; % Default number of means to use to abinitio reconstruction.
message=sprintf('Number of class means to use to abinitio reconstruction of each group? ');
nmeans=fmtinput(message,defnmeans,'%d');


%% Update workflow struct

workflow.abinitio.nmeans=nmeans; % Number of means to use to abinitio reconstruction.

tree=struct2xml(workflow);
save(tree,workflow_fname); 

%% Generate abinitio models
cryo_workflow_abinitio_exectue(workflow_fname);
