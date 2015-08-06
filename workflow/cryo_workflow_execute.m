function cryo_workflow_execute(workflow_fname)
% CRYO_WORKFLOW_EXECUTE     Execute complete reconstruction workflow
%
% cryo_workflow_execute(workflow_fname)
%   Run all reconstruction steps with parameters given in the file
%   workflow_fname. The steps are preprocessing, classification, generating
%   class means, and abinitio reconstruction.
%
% Yoel Shkolnisky, August 2015.

cryo_workflow_preprocess_execute(workflow_fname);
cryo_workflow_classify_execute(workflow_fname);
cryo_workflow_classmeans_execute(workflow_fname);
cryo_workflow_abinitio_execute(workflow_fname);


