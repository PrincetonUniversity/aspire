function cryo_workflow_classmeans_execute(workflow_fname)
% CRYO_WORKFLOW_CLASSMEANS_EXECUTE  Generate class means
%
% cryo_workflow_classmeans_execute(workflow_fname)
%   Generate class means using precomputed classification according to the
%   parameters stored in the file workflow_fname.
%
% See also cryo_workflow_classmeans
%
% Yoel Shkolnisky, August 2015.

%% Validate workflow file
cryo_workflow_classmeans_validate(workflow_fname);

%% Read workflow file
tree=xmltree(workflow_fname);
workflow=convert(tree);

%% Execute preprocessing

open_log(fullfile(workflow.info.working_dir,workflow.info.logfile));

nnavg=str2double(workflow.classmeans.nnavg);
numgroups=str2double(workflow.preprocess.numgroups); 

for groupid=1:numgroups   
    fname=sprintf('phaseflipped_cropped_downsampled_prewhitened_group%d.mrc',groupid);
    fullfilename=fullfile(workflow.info.working_dir,fname);
    log_message('Loading prewhitened projection from %s',fname);
    prewhitened_projs=ReadMRC(fullfilename);
    
    matname=fullfile(workflow.info.working_dir,sprintf('VDM_data_%d',groupid));
    log_message('Loading %s',matname);
    load(matname);
    
    log_message('Starting align_main');
    list_recon=1:size(prewhitened_projs,3);
    [ shifts, corr, average, norm_variance ] = align_main(prewhitened_projs,...
        VDM_angles, class_VDM, class_VDM_refl, FBsPCA_data,...
            nnavg,15, list_recon);
    [~,classcoreidx]=sort(norm_variance); % classcoreidx are the
            % indices of the most consistent class averages. The
            % corresponding phaseflipped images will be used for
             % reconstruction.
    log_message('Finished align_main');
    
    % Save averages sorted by norm variance
    average=average(:,:,classcoreidx);
    [average,doflip]=cryo_globalphaseflip(average); % Check if a global phase flip should be applied
    
    if doflip
        log_message('Applying global phase flip to averages');
    end
    
    fname=sprintf('averages_nn%02d_group%d.mrc',nnavg,groupid);
    fname=fullfile(workflow.info.working_dir,fname);
    log_message('Saving %s',fname);
    WriteMRC(single(average),1,fname);
    
    reloadname=sprintf('averages_info_nn%02d_group%d',nnavg,groupid);
    save(fullfile(workflow.info.working_dir,reloadname),...
        'shifts','corr','norm_variance','classcoreidx');
    
    % Compute effetive CTF of each average    
    fname=sprintf('ctfs_group%d.star',groupid);
    fullfilename=fullfile(workflow.info.working_dir,fname);    
    log_message('Load %s',fullfilename);
    CTFdata=readSTAR(fullfilename);
    n=size(average,1);
    effectiveCTFs=zeros(size(average));
    
    log_message('Computing effective CTFs for group %d',groupid);
    printProgressBarHeader;
    
    for k=1:size(average,3)
        progressTicFor(k,size(average,3));
        idx=classcoreidx(k); % Index of the average in unsorted stack of averages
        ectf=zeros(n);      % Effective CTF for the current average.
        for nnk=1:nnavg
            nnidx=class_VDM(idx,nnk);
            [voltage,DefocusU,DefocusV,DefocusAngle,Cs,pixA,A]=...
                cryo_parse_Relion_CTF_struct(CTFdata.data{nnidx});
            h=cryo_CTF_Relion(n,voltage,DefocusU,DefocusV,DefocusAngle,...
                Cs,pixA,A);
            if str2double(workflow.preprocess.phaseflip)
                h=abs(h);
            end
            ectf=ectf+h;
        end
        effectiveCTFs(:,:,k)=ectf./nnavg;
    end
    
    log_message('Saving effective CTFs for group %d',groupid);
    fname=sprintf('ctfs_effective_nn%02d_group%d.mrc',nnavg,groupid);
    fullfilename=fullfile(workflow.info.working_dir,fname);
    WriteMRC(single(effectiveCTFs),1,fullfilename);
end

log_message('Workflow file: %s\n',workflow_fname);
log_message('Use this file name when calling subsequent funtions.\n');
log_message('Call next cryo_workflow_abinitio(''%s'')\n',workflow_fname);

