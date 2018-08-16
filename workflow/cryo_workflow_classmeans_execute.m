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

log_message('Starting cryo_workflow_classmeans_execute');
log_message('Loaded XML file %s (MD5: %s)',workflow_fname,MD5(workflow_fname));

nnavg=str2double(workflow.classmeans.nnavg);
numgroups=str2double(workflow.preprocess.numgroups);

use_EM=str2double(workflow.classmeans.use_EM);
gpu_list=-1;
if use_EM
    gpu_list=str2double(workflow.classmeans.gpu_list);
    num_EM_averages=str2double(workflow.classmeans.num_EM_averages);
end


for groupid=1:numgroups
    fname=sprintf('phaseflipped_cropped_downsampled_prewhitened_group%d.mrcs',groupid);
    fullfilename=fullfile(workflow.info.working_dir,fname);
    
    log_message('Loading prewhitened projections from %s (MD5: %s)',fname,MD5(fullfilename));

    % Since we are loading all nearest neighbors of an image, there will be
    % many cache misses. So no point in have a large chache.
    prewhitened_projs=imagestackReader(fullfilename,1);
    
    matname=fullfile(workflow.info.working_dir,sprintf('VDM_data_%d.mat',groupid));
    log_message('Loading %s (MD5: %s)',matname,MD5(matname));
    classification_data=load(matname);
    
    log_message('Starting align_main');
    [~,tmpkey]=fileparts(tempname); % Generate a unique temporary directory name
    tmpdir=fullfile(tempmrcdir,'align_main',tmpkey);
    log_message('Using temporary directory %s',tmpdir);
    if ~exist(tmpdir,'dir')
        mkdir(tmpdir);
    end
    delete(fullfile(tmpdir,'*')); % Delete any leftovers from the temp directory

    list_recon=1:prewhitened_projs.dim(3); % Compute class averages for all input images
    log_message('Generating %d class averages. use_EM=%d',numel(list_recon),use_EM);
    [ shifts, corr, unsortedaveragesfname, ~, norm_variance] = align_main( prewhitened_projs,...
        classification_data.VDM_angles, classification_data.class_VDM,...
        classification_data.class_VDM_refl, classification_data.sPCA_data,...
        nnavg, 15, list_recon, 0, [],tmpdir);
    log_message('Finished align_main');         

    % Write unsorted stacks
    fnameunsorted=sprintf('averages_nn%02d_unsorted_group%d.mrcs',nnavg,groupid);
    fnameunsorted=fullfile(workflow.info.working_dir,fnameunsorted);
    
    % Save unsorted averages and apply global phaseglip if needed.
    doflip=cryo_globalphaseflip_outofcore(unsortedaveragesfname ,fnameunsorted);
    if doflip
        log_message('Global contrast inversion applied to unsorted averages.');
    end

    % Sort the resulting class averages by their contrast.
    averages_contrast=cryo_image_contrast_outofcore(fnameunsorted);
    [~,classcoreidx]=sort(averages_contrast,'descend'); %classcoreidx are 
        % the indices of the class averages sorted from the most consistent
        % average to the least consistent one.
    fnamesorted=sprintf('averages_nn%02d_sorted_group%d.mrcs',nnavg,groupid);
    fnamesorted=fullfile(workflow.info.working_dir,fnamesorted);
    log_message('Writing sorted averages (by contrast) into %s',fnamesorted);
    cryo_sort_stack_outofcore(fnameunsorted,classcoreidx,fnamesorted);
    
    % Print how many unique raw images are involved in the top i'th
    % percentile of class averages (sorted by averages_contrast). Also
    % print for each row the norm_variance of the i'th percentile class
    % average. 
    log_message('Number of raw projection in each the each percentile of averages (percentile,num means, num unique raw, contrast)');
    nprojs=numel(list_recon);
    for k=1:10
        nmeans=round(nprojs*k/10);
        ii=classification_data.class_VDM(classcoreidx(1:nmeans),:);
        nrawprojs=numel(unique(ii));
        log_message('\t%3d%%\t %7d \t %7d \t %4.2e',k*10,nmeans,nrawprojs,averages_contrast(classcoreidx(nmeans)));
    end
        
    
    if use_EM
    % To avoid major rewrites at this time, I am just calling align_main
    % twice. Once for all images without EM and once for a selected subset
    % using EM. It is cleaner to separate align_main from the EM step, but
    % I wanted to make the change incrementaly to make the change more
    % quickly. Yoel Shkolnisky August 2018
    
    % Select a subset of the class averages to refine using EM.
    % We are working with the sorted class averages since they are sorted by
    % quality. Since the classification data is given for the unsorted class
    % averages, we need to changes the indices in the nearest neighbors table
    % to the sorted order. The mapping from the old order to the new order is
    % enocded by classcoreidx. The image corresponding to image k in the sorted
    % order is classcoreidx(k).
    
    NNtbl=classification_data.class_VDM(classcoreidx,:);
    NNtblreordered=NNtbl;
    % Fix the indices in NNtbl to correspond to the sorted order of the images.
    imap=zeros(1,max(classcoreidx));
    imap(classcoreidx)=1:numel(classcoreidx);
    for k=1:numel(NNtbl)
        NNtblreordered(k)=imap(NNtbl(k));
    end
    
    % Select a subset of the averages to refine
    list_EM_recon=cryo_select_subset(NNtblreordered,num_EM_averages,min(size(NNtblreordered,1),20000));
    % Select a subset only from the first 2000 images. This is a temporary
    % hack. The user should choose that according to the quality of the
    % class averages.
    list_EM_recon=classcoreidx(list_EM_recon); % Map back to indices in the unsorted file of averages.
    
    % Save an MRCS file with the averages to refine
    em_seeds_fname=sprintf('em_seeds_nn%02d_group%d.mrcs',nnavg,groupid);
    unsortedstackreader=imagestackReader(fnameunsorted);
    outstack=imagestackWriter(em_seeds_fname,num_EM_averages);
    for k=1:num_EM_averages
        proj=unsortedstackreader.getImage(list_EM_recon(k));
        outstack.append(proj);
    end
    outstack.close;
    
    delete(fullfile(tmpdir,'*')); % Delete leftovers from previous alignment
    [ ~, ~, ~, unsortedaveragesEMfname,~,loglikelihood] = align_main( prewhitened_projs,...
        classification_data.VDM_angles, classification_data.class_VDM,...
        classification_data.class_VDM_refl, classification_data.sPCA_data,...
        nnavg, 15, list_EM_recon, use_EM, gpu_list,tmpdir);
 
        fnameunsorted=sprintf('averages_nn%02d_EM_group%d.mrcs',nnavg,groupid);
        fnameunsorted=fullfile(workflow.info.working_dir,fnameunsorted);
        doflip=cryo_globalphaseflip_outofcore(unsortedaveragesEMfname ,fnameunsorted);
        if doflip
            log_message('Global contrast inversion applied to unsorted EM averages.');
        end
        
    end
    
    reloadname=sprintf('averages_info_nn%02d_group%d.mat',nnavg,groupid);
    save(fullfile(workflow.info.working_dir,reloadname),'-v7.3',...
        'shifts','corr','averages_contrast','classcoreidx',...
        'classification_data','doflip', 'loglikelihood');
    
    log_message('Saving %s (MD5: %s)',fullfile(workflow.info.working_dir,reloadname)...
        ,MD5(fullfile(workflow.info.working_dir,reloadname)));
    
    delete(fullfile(tmpdir,'*')); % Clean the temp directory.
    %
    %     % Compute effetive CTF of each average
    %     fname=sprintf('ctfs_group%d.star',groupid);
    %     fullfilename=fullfile(workflow.info.working_dir,fname);
    %     if exist(fullfilename,'file')~=2
    %         log_message('%s does not exist. Skipping computing effective CTF of each class average',fullfilename);
    %     else
    %
    %         log_message('Load %s',fullfilename);
    %         CTFdata=readSTAR(fullfilename);
    %         n=size(average,1);
    %         effectiveCTFs=zeros(size(average));
    %
    %         log_message('Computing effective CTFs for group %d',groupid);
    %         printProgressBarHeader;
    %
    %         for k=1:size(average,3)
    %             progressTicFor(k,size(average,3));
    %             idx=classcoreidx(k); % Index of the average in unsorted stack of averages
    %             ectf=zeros(n);      % Effective CTF for the current average.
    %             for nnk=1:nnavg
    %                 nnidx=class_VDM(idx,nnk);
    %                 [voltage,DefocusU,DefocusV,DefocusAngle,Cs,pixA,A]=...
    %                     cryo_parse_Relion_CTF_struct(CTFdata.data{nnidx});
    %                 h=cryo_CTF_Relion(n,voltage,DefocusU,DefocusV,DefocusAngle,...
    %                     Cs,pixA,A);
    %                 if str2double(workflow.preprocess.phaseflip)
    %                     h=abs(h);
    %                 end
    %                 ectf=ectf+h;
    %             end
    %             effectiveCTFs(:,:,k)=ectf./nnavg;
    %         end
    %
    %         log_message('Saving effective CTFs for group %d',groupid);
    %         fname=sprintf('ctfs_effective_nn%02d_group%d.mrcs',nnavg,groupid);
    %         fullfilename=fullfile(workflow.info.working_dir,fname);
    %         WriteMRC(single(effectiveCTFs),1,fullfilename);
    %     end
end

log_message('Workflow file: %s\n',workflow_fname);
log_message('Use this file name when calling subsequent funtions.\n');
log_message('Call next cryo_workflow_abinitio(''%s'')\n',workflow_fname);

close_log;
