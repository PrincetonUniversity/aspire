function cryo_workflow_classify(workflow_fname)

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