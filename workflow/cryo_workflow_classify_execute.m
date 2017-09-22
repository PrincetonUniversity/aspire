function cryo_workflow_classify_execute(workflow_fname)
% CRYO_WORKFLOW_CLASSIFY_EXECUTE  Execute data set classification
%
% cryo_workflow_classify_execute(workflow_fname)
%   Classify preprocessed projections according to the parameters stored
%   in the file workflow_fname.
%
% See also cryo_workflow_classify
%
% Yoel Shkolnisky, August 2015.

%% Validate workflow file
cryo_workflow_classify_validate(workflow_fname);

%% Read workflow file
tree=xmltree(workflow_fname);
workflow=convert(tree);

%% Execute classification
open_log(fullfile(workflow.info.working_dir,workflow.info.logfile));

log_message('Starting cryo_workflow_classify_execute');
log_message('Loaded XML file %s (MD5: %s)',workflow_fname,MD5(workflow_fname));

numgroups=str2double(workflow.preprocess.numgroups); 

for groupid=1:numgroups
    % Read prewhitened projections
    fname=sprintf('phaseflipped_cropped_downsampled_prewhitened_group%d.mrc',groupid);
        
    fullfilename=fullfile(workflow.info.working_dir,fname);
    log_message('Loading %s (MD5: %s)',fullfilename,MD5(fullfilename));
    prewhitened_projs=ReadMRC(fullfilename);
    n=size(prewhitened_projs,1);
    prewhitened_projs=double(prewhitened_projs); % Convert to double for VDM below.
    
    % Estimate variance of the noise. Should be 1 (if images have been normalized).
    log_message('Estimating SNR of images');
    [~,~,var_n]=cryo_estimate_snr(prewhitened_projs);
    log_message('Estimated noise variance=%d',var_n);    
    
    
    log_message('Starting computing steerable PCA');
    sPCA_data=data_sPCA(prewhitened_projs,  var_n);
    log_message('Finished computing steerable PCA');
    
    log_message('Starting class averaging initial classificaiton');
    [ class, class_refl, rot, ~,  ~] = Initial_classification_FD_update(sPCA_data,...
        str2double(workflow.classification.n_nbor),...
        str2double(workflow.classification.isrann));    
    log_message('Finished class averaging initial classificaiton');
    
    log_message('Starting VDM');
    [ class_VDM, class_VDM_refl, VDM_angles ] = VDM(class, ones(size(class)), rot,...
        class_refl, str2double(workflow.classification.k_VDM_in),...
        str2double(workflow.classification.VDM_flag),...
        str2double(workflow.classification.k_VDM_out));
    log_message('Finished VDM classification...');
    
    matname=fullfile(workflow.info.working_dir,sprintf('VDM_data_%d.mat',groupid));
    save(matname,'class','class_refl','rot','sPCA_data','class_VDM',...
        'class_VDM_refl','class_VDM_refl','VDM_angles');
    
    log_message('Saved %s (MD5: %s)',matname,MD5(matname));
end

log_message('Workflow file: %s',workflow_fname);
log_message('Use this file name when calling subsequent funtions.');
log_message('Call next cryo_workflow_classmeans(''%s'')',workflow_fname);

close_log;
