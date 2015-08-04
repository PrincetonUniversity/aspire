function cryo_workflow_classify_execute(workflow_fname)

%% Validate workflow file
cryo_workflow_classify_validate(workflow_fname);

%% Read workflow file
tree=xmltree(workflow_fname);
workflow=convert(tree);

%% Execute classification
open_log(fullfile(workflow.info.working_dir,workflow.info.logfile));

numgroups=str2double(workflow.preprocess.numgroups); 

for groupid=1:numgroups
    % Read prewhitened projections
    fname=sprintf('phaseflipped_downsampled_prewhitened_group%d.mrc',groupid);
    fullfilename=fullfile(workflow.info.working_dir,fname);
    prewhitened_projs=ReadMRC(fullfilename);
    n=size(prewhitened_projs,1);
    prewhitened_projs=double(prewhitened_projs); % Convert to double for VDM below.
    
    log_message('Starting class averaging initial classificaiton');
    r_max=round(n/2)-10;
    [ class, class_refl, rot, ~, FBsPCA_data, ~ ] = Initial_classification(prewhitened_projs , r_max,...
        str2double(workflow.classification.n_nbor),...
        str2double(workflow.classification.isrann));
    log_message('Finished class averaging initial classificaiton');
    
    log_message('Starting VDM');
    [ class_VDM, class_VDM_refl, VDM_angles ] = VDM(class, ones(size(class)), rot,...
        class_refl, str2double(workflow.classification.k_VDM_in),...
        str2double(workflow.classification.VDM_flag),...
        str2double(workflow.classification.k_VDM_out));
    disp('Finished VDM classification...');
    
    matname=fullfile(workflow.info.working_dir,sprintf('VDM_data_%d',groupid));
    save(matname,'class','class_refl','rot','FBsPCA_data','class_VDM',...
        'class_VDM_refl','class_VDM_refl','VDM_angles');
end

log_message('Workflow file: %s\n',workflow_fname);
log_message('Use this file name when calling subsequent funtions.\n');
log_message('Call next cryo_workflow_classmeans(''%s'')\n',workflow_fname);

close_log;