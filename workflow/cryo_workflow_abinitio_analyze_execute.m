function cryo_workflow_abinitio_analyze_execute(workflow_fname)
% CRYO_WORKFLOW_ANALYZE_EXECUTE  Execute analysis of reconstructed abinitio models
%
% cryo_workflow_abinitio_analyze_execute(workflow_fname)
%   Analyze the abinitio reconstructed models. The analysis includes
%   aligning the volumes, plotting FSC curves, and plotting viewing
%   directions.
%
% See also cryo_workflow_abinitio_analyze
%
% Yoel Shkolnisky, October 2015.

% Validate that the given workflow has all required paramters.
cryo_workflow_abinitio_analyze_validate(workflow_fname);

%% Read workflow file
tree=xmltree(workflow_fname);
workflow=convert(tree);

% Set workflow parameters
numgroups=str2double(workflow.preprocess.numgroups);
nmeans=str2double(workflow.abinitio.nmeans);
nnavg=str2double(workflow.abinitio.nnavg);
pixA=str2double(workflow.analysis.pixA);

%% Align volumes reconstructed from averages
log_message('Align volumes reconstructed from averages');

% Load first volume
vol1fname=sprintf('vol_nn%d_nm%d_group%d.mrc',nnavg,nmeans,1);
vol1fname=fullfile(workflow.info.working_dir,vol1fname);
log_message('Loading %s',vol1fname);
vol1=ReadMRC(vol1fname);

% Align all remainging volumes agains the first.
for groupid=2:numgroups
    % Load second volume
    vol2fname=sprintf('vol_nn%d_nm%d_group%d.mrc',nnavg,nmeans,groupid);
    vol2fname=fullfile(workflow.info.working_dir,vol2fname);
    log_message('Loading %s',vol2fname);    
    vol2=ReadMRC(vol2fname);
    
    % Align vol1 and vol2 
    [~,~,vol2aligned]=cryo_align_densities(vol1,vol2,pixA,1);

    % Save aligned volume
    vol2fname=sprintf('vol_nn%d_nm%d_group%d_aligned.mrc',nnavg,nmeans,groupid);
    vol2fname=fullfile(workflow.info.working_dir,vol2fname);
    WriteMRC(vol2aligned,1,vol2fname);

    % Plot FSC
    [resA,h]=plotFSC(vol1,vol2aligned,0.143,pixA);
    log_message('FSC resolution of group 1 and group %d is %d Angstroms',...
        groupid,resA);
    
    fscFIGname=fullfile(workflow.info.working_dir,...
        sprintf('fsc_nn%d_nm_%d_group1to%d.fig',nnavg,nmeans,groupid));
    fscEPSname=fullfile(workflow.info.working_dir,...
        sprintf('fsc_nn%d_nm_%d_group1to%d.eps',nnavg,nmeans,groupid));
    hgsave(fscFIGname);
    print('-depsc',fscEPSname);
    close(h);
    
end

%% Align volumes reconstructed from phaseflipped raw projections
log_message('Align volumes reconstructed from raw projections');

% Load first volume
vol1fname=sprintf('vol_raw_nn%d_nm%d_group%d.mrc',nnavg,nmeans,1);
vol1fname=fullfile(workflow.info.working_dir,vol1fname);
log_message('Loading %s',vol1fname);
vol1=ReadMRC(vol1fname);

% Align all remainging volumes agains the first.
for groupid=2:numgroups
    % Load second volume
    vol2fname=sprintf('vol_raw_nn%d_nm%d_group%d.mrc',nnavg,nmeans,groupid);   
    vol2fname=fullfile(workflow.info.working_dir,vol2fname);
    log_message('Loading %s',vol2fname);
    vol2=ReadMRC(vol2fname);
    
    % Align vol1 and vol2 
    [~,~,vol2aligned]=cryo_align_densities(vol1,vol2,pixA,1);

    % Save aligned volume
    vol2fname=sprintf('vol_raw_nn%d_nm%d_group%d_aligned.mrc',nnavg,nmeans,groupid);
    vol2fname=fullfile(workflow.info.working_dir,vol2fname);
    WriteMRC(vol2aligned,1,vol2fname);

    % Plot FSC
    [resA,h]=plotFSC(vol1,vol2aligned,0.143,pixA);
    log_message('FSC resolution of group 1 and group %d is %d Angstroms',...
        groupid,resA);
    
    fscFIGname=fullfile(workflow.info.working_dir,...
        sprintf('fsc_raw_nn%d_nm_%d_group1to%d.fig',nnavg,nmeans,groupid));
    fscEPSname=fullfile(workflow.info.working_dir,...
        sprintf('fsc_raw_nn%d_nm_%d_group1to%d.eps',nnavg,nmeans,groupid));
    hgsave(fscFIGname);
    print('-depsc',fscEPSname);
    close(h);
    
end

% %% Align CTF corrected volumes
% 
% % Load first volume
% vol1fname=sprintf('vol_ctf_corrected_nn%d_nm%d_group%d.mrc',nnavg,nmeans,1);
% vol1fname=fullfile(workflow.info.working_dir,vol1fname);
% 
% if exist(vol1fname,'file')~=2
%     log_message('%s does not exist. Skipping analyzing CTF corrected volumes',vol1fname);
% else
%     
%     log_message('Align CTF corrected volumes');
%     vol1=ReadMRC(vol1fname);
%     
%     % Align all remainging volumes agains the first.
%     for groupid=2:numgroups
%         % Load second volume
%         vol2fname=sprintf('vol_ctf_corrected_nn%d_nm%d_group%d.mrc',nnavg,nmeans,groupid);
%         vol2fname=fullfile(workflow.info.working_dir,vol2fname);
%         vol2=ReadMRC(vol2fname);
%         
%         % Align vol1 and vol2
%         [~,~,vol2aligned]=cryo_align_densities(vol1,vol2,pixA,1);
%         
%         % Save aligned volume
%         vol2fname=sprintf('vol_ctf_corrected_nn%d_nm%d_group%d_aligned.mrc',nnavg,nmeans,groupid);
%         vol2fname=fullfile(workflow.info.working_dir,vol2fname);
%         WriteMRC(vol2aligned,1,vol2fname);
%         
%         % Plot FSC
%         [resA,h]=plotFSC(vol1,vol2aligned,0.143,pixA);
%         log_message('FSC resolution of group 1 and group %d is %d Angstroms',...
%             groupid,resA);
%         
%         fscFIGname=fullfile(workflow.info.working_dir,...
%             sprintf('fsc_ctf_corrected_nn%d_nm_%d_group1to%d.fig',nnavg,nmeans,groupid));
%         fscEPSname=fullfile(workflow.info.working_dir,...
%             sprintf('fsc_ctf_corrected_nn%d_nm_%d_group1to%d.eps',nnavg,nmeans,groupid));
%         hgsave(fscFIGname);
%         print('-depsc',fscEPSname);
%         close(h);        
%     end
% end

%% Plot viewing directions
log_message('Plot viewing directions');

for groupid=2:numgroups
    matname=sprintf('abinitio_info_nn%d_nm%d_group%d',nnavg,nmeans,groupid);
    matname=fullfile(workflow.info.working_dir,matname);
    s=load(matname,'rotations');
    [~,h2]=cryo_plot_viewing_directions(s.rotations);
    figure(h2);
    rotFIGname=fullfile(workflow.info.working_dir,...
        sprintf('rotations_nn%d_nm_%d_group1to%d.fig',nnavg,nmeans,groupid));
    rotEPSname=fullfile(workflow.info.working_dir,...
        sprintf('rotations_nn%d_nm_%d_group1to%d.eps',nnavg,nmeans,groupid));
    hgsave(rotFIGname);
    print('-depsc',rotEPSname);
    close(h2);
end

