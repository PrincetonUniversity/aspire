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

rand('twister', 1137);
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
    % Clear variables from previous groups   
    clear sPCA_data class class_refl rot class_VDM class_VDM_refl VDM_angles
    
    % Read prewhitened projections
    fname=sprintf('phaseflipped_cropped_downsampled_prewhitened_group%d.mrcs',groupid);
        
    fullfilename=fullfile(workflow.info.working_dir,fname);
    log_message('Loading %s (MD5: %s)',fullfilename,MD5(fullfilename));
    prewhitened_projs=ReadMRC(fullfilename);
    prewhitened_projs=double(prewhitened_projs); % Convert to double for VDM below.
    
    % Estimate variance of the noise. Should be 1 (if images have been normalized).
    log_message('Estimating SNR of images');
    [snr,var_s,var_n]=cryo_estimate_snr(prewhitened_projs);
    log_message('Estimated SNR=1/%d (more precisely %d)',round(1/snr),snr);    
    log_message('Estimated signal variance=%d',var_s);    
    log_message('Estimated noise variance=%d',var_n);    
        
    log_message('Starting computing steerable PCA');
    % sPCA_data_fb=data_sPCA(prewhitened_projs,  var_n);
    % log_message('PSWF sPCA');
    % sPCA_data_pswf=sPCA_PSWF(prewhitened_projs,  var_n);
    sPCA_data=sPCA_PSWF(prewhitened_projs,  var_n);
    log_message('Finished computing steerable PCA');
    
    % sPCA_data=sPCA_data_fb;
    % sPCA_data=sPCA_data_pswf;
    
    % Janes original code did not use all sPCA components for initial
    % classification. To save memory, it used up to 400 components. I
    % experimented with 100, 200, 400, and 800 components, and it seems
    % that using fewer components results in better class averages, in the
    % sense that there are more class averages of good quality. However,
    % from my limited tests, it seems that the effect on the final
    % reconstruction (if one uses say 200 highest-contrast averages) is
    % negligible. I therefore use some default value for ncomp in
    % cryo_workflow_classify, and let the user change it. 
    % Yoel Shkolnisky, July 2018.
    ncomp=str2double(workflow.classification.ncomp); 
    ncomp=min(ncomp,numel(sPCA_data.eigval));
    log_message('sPCA resulted in %d components',numel(sPCA_data.eigval));
    log_message('Using %d sPCA components for classification',ncomp);
    
    sPCA_data.eigval=sPCA_data.eigval(1:ncomp);
    sPCA_data.Freqs=sPCA_data.Freqs(1:ncomp);
    sPCA_data.RadFreqs=sPCA_data.RadFreqs(1:ncomp);
    sPCA_data.Coeff=sPCA_data.Coeff(1:ncomp,:);
    sPCA_data.eig_im=sPCA_data.eig_im(:,1:ncomp);
    matname=fullfile(workflow.info.working_dir,sprintf('VDM_data_%d.mat',groupid));
    save(matname,'-v7.3','sPCA_data');

    
    log_message('Starting class averaging initial classificaiton');
    [ class, class_refl, rot, ~,  ~] = Initial_classification_FD_update(sPCA_data,...
        str2double(workflow.classification.n_nbor),...
        str2double(workflow.classification.isrann));    
    log_message('Finished class averaging initial classificaiton');

    save(matname,'class','class_refl','rot','-append');    
    
    log_message('Starting VDM');
    [ class_VDM, class_VDM_refl, VDM_angles ] = VDM(class, ones(size(class)), rot,...
        class_refl, str2double(workflow.classification.k_VDM_in),...
        str2double(workflow.classification.VDM_flag),...
        str2double(workflow.classification.k_VDM_out));
    log_message('Finished VDM classification...');

save(matname,'class_VDM','class_VDM_refl','VDM_angles','-append');
    
    log_message('Saved %s (MD5: %s)',matname,MD5(matname));
end

log_message('Workflow file: %s',workflow_fname);
log_message('Use this file name when calling subsequent funtions.');
log_message('Call next cryo_workflow_classmeans(''%s'')',workflow_fname);

close_log;
