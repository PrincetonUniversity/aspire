function cryo_workflow_classify(workflow_fname)
% CRYO_WORKFLOW_CLASSIFY  Interactive data set classification
%
% cryo_workflow_classify(workflow_fname)
%   Interactively collect all parameters required to classify preprocessed
%   projections data set and execute the classification. No class means are
%   produced at this step. Call cryo_workflow_classmeans to generate class
%   means using the computed classification.
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


%% Update workflow struct

workflow.classification.n_nbor=100; %number of nearest neighbors for initial classification.
workflow.classification.isrann = 0;
workflow.classification.k_VDM_in = 20; % number of nearest neighbors for building graph for VDM.
workflow.classification.VDM_flag = 0;
workflow.classification.k_VDM_out = 100; % output number of nearest neighbors

tree=struct2xml(workflow);
save(tree,workflow_fname); 


%% Execute classification
cryo_workflow_classify_execute(workflow_fname);