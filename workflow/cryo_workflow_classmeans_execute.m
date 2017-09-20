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

for groupid=1:numgroups
    fname=sprintf('phaseflipped_cropped_downsampled_prewhitened_group%d.mrc',groupid);
    fullfilename=fullfile(workflow.info.working_dir,fname);
    
    log_message('Loading prewhitened projections from %s (MD5: %s)',fname,MD5(fullfilename));

    %prewhitened_projs=ReadMRC(fullfilename);
    
    % Since we are loading all nearest neighbors of an image, there will be
    % many cache misses. So no point in have a large chache.
    prewhitened_projs=imagestackReader(fullfilename,1);
    
    matname=fullfile(workflow.info.working_dir,sprintf('VDM_data_%d.mat',groupid));
    log_message('Loading %s (MD5: %s)',matname,MD5(matname));
    load(matname);
    
    log_message('Starting align_main');
    tmpdir=fullfile(workflow.info.working_dir,'tmp');
    if ~exist(tmpdir,'dir')
        mkdir(tmpdir);
    end
    delete(fullfile(tmpdir,'*')); % Delete any leftovers from the temp directory

    list_recon=1:prewhitened_projs.dim(3);
    [ shifts, corr, unsortedaveragesfname, norm_variance ] = align_main( prewhitened_projs,...
        VDM_angles, class_VDM, class_VDM_refl, sPCA_data, nnavg, 15, list_recon, tmpdir);
    log_message('Finished align_main');         

    
%     [~,classcoreidx]=sort(norm_variance); % classcoreidx are the
%     % indices of the most consistent class averages. The
%     % corresponding phaseflipped images will be used for
%     % reconstruction.

    % Sort the resulting class averages by their contrast.
    averages_contrast=cryo_image_contrast_outofcore(unsortedaveragesfname);
    [~,classcoreidx]=sort(averages_contrast,'descend'); % Only averages with highest 
    %  will be used for reconstruction.
    
    
    % Print how many unique raw images are involved in the top i'th
    % percentile of class averages (sorted by averages_contrast). Also
    % print for each row the norm_variance of the i'th percentile class
    % average. 
    log_message('Number of raw projection in each the each percentile of averages (percentile,num means, num unique raw, contrast)');
    nprojs=prewhitened_projs.dim(3);
    for k=1:10
        nmeans=round(nprojs*k/10);
        ii=class_VDM(classcoreidx(1:nmeans),:);
        nrawprojs=numel(unique(ii));
        log_message('\t%3d%%\t %7d \t %7d \t %4.2e',k*10,nmeans,nrawprojs,averages_contrast(classcoreidx(nmeans)));
    end
    
    % Determine if global phase flip is required
    log_message('Checking if global phase flip is required.');
    [doflip,signalmean,noisemean]=cryo_globalphaseflip_outofcore(unsortedaveragesfname); % Check if a global phase flip should be applied    
    log_message('Signal mean = %5.3e, noise power = %5.3e',signalmean,noisemean);
    if doflip
        log_message('Applied global phase flip to averages');
    end
    flipflag=(-2)*doflip+1;
    
    % Save averages sorted by norm variance    
    fname=sprintf('averages_nn%02d_group%d.mrc',nnavg,groupid);
    fname=fullfile(workflow.info.working_dir,fname);
    log_message('Writing sorted averages (by contrast) into %s',fname);
            
    average=imagestackReader(unsortedaveragesfname,100);
    sortedaverages=imagestackWriter(fname,average.dim(3),1,100);
    for k=1:numel(classcoreidx)
        im=average.getImage(classcoreidx(k));
        im=im.*flipflag;
        sortedaverages.append(im);        
    end
    sortedaverages.close;
    
    log_message('Saved averages in %s (MD5: %s)',fname,MD5(fname))
%    average=average(:,:,classcoreidx);
%    [average,doflip]=cryo_globalphaseflip(average); % Check if a global phase flip should be applied
   
%     fname=sprintf('averages_nn%02d_group%d.mrc',nnavg,groupid);
%     fname=fullfile(workflow.info.working_dir,fname);
%     log_message('Saving %s',fname);
%     WriteMRC(single(average),1,fname);
    
    reloadname=sprintf('averages_info_nn%02d_group%d.mat',nnavg,groupid);
    save(fullfile(workflow.info.working_dir,reloadname),...
        'shifts','corr','averages_contrast','classcoreidx','VDM_angles',...
        'class_VDM', 'class_VDM_refl','doflip');
    
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
%         fname=sprintf('ctfs_effective_nn%02d_group%d.mrc',nnavg,groupid);
%         fullfilename=fullfile(workflow.info.working_dir,fname);
%         WriteMRC(single(effectiveCTFs),1,fullfilename);
%     end
end

log_message('Workflow file: %s\n',workflow_fname);
log_message('Use this file name when calling subsequent funtions.\n');
log_message('Call next cryo_workflow_abinitio(''%s'')\n',workflow_fname);

close_log;
