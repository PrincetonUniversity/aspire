function cryo_workflow_abinitio_analyze(workflow_fname)
% CRYO_WORKFLOW_ABINITIO_ANALYZE  Analyze the reconstructed abinitio models
%
% cryo_workflow_abinitio_analyze(workflow_fname)
%   Analyze the abinitio models reconstructed by the workflow given by
%   workflow_name. The analysis consists of aligning and saving the
%   reconstructed abinitio models, plotting FSC curves, and ploting the
%   estiated viewing directins.
%
% Yoel shkolnisky, September 2015.


%% Read workflow file
tree=xmltree(workflow_fname);
workflow=convert(tree);

%% Get required parameters.
% If pixAdownsampled exists, then use it.
% Otherwise, read it from the user.

if isfield(workflow.preprocess,'pixAdownsampled') 
    pixA=workflow.preprocess.pixAdownsampled;
else
    message='Pixel size (Angstrom)?';
    pixA=fmtinput(message,'','%d');
end

%% Update workflow struct
workflow.analysis.pixA=pixA;

tree=struct2xml(workflow);
save(tree,workflow_fname); 

%% Execute analysis
cryo_workflow_abinitio_analyze_execute(workflow_fname);
