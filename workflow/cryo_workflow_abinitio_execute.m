function cryo_workflow_abinitio_execute(workflow_fname)
% CRYO_WORKFLOW_CLASSIFY_EXECUTE  Reconstruct abinitio models
%
% cryo_workflow_abinitio_execute(workflow_fname)
%   Generate abinitio models from precomputed class averages according to
%   the parameters stored in the file workflow_fname.
%
% See also cryo_workflow_abinitio
%
% Yoel Shkolnisky, August 2015.

%% Validate workflow file
cryo_workflow_abinitio_validate(workflow_fname);

%% Read workflow file
tree=xmltree(workflow_fname);
workflow=convert(tree);

%% Reconstruct abinitio models

open_log(fullfile(workflow.info.working_dir,workflow.info.logfile));

numgroups=str2double(workflow.preprocess.numgroups);
nmeans=str2double(workflow.abinitio.nmeans);
nnavg=str2double(workflow.abinitio.nnavg);

for groupid=1:numgroups    
    reloadname=sprintf('averages_info_nn%02d_group%d',nnavg,groupid);
    reloadname=fullfile(workflow.info.working_dir,reloadname);
    log_message('Loading %s',reloadname);
    load(reloadname);
    fname=sprintf('averages_nn%02d_group%d.mrc',nnavg,groupid);
    average=ReadMRC(fullfile(workflow.info.working_dir,fname));
    average=average(:,:,1:nmeans);
    
    K=size(average,3);
    log_message('Averages loaded. Using K=%d averages of size %d x %d',K,size(average,1),size(average,2));
    
    % Mask projections    
    mask_radius=round(size(average,1)*0.45);
    log_message('Masking radius is %d pixels',mask_radius);
    [masked_average,~]=mask_fuzzy(average,mask_radius);

    log_message('Computeing polar Fourier transform of class averages');
    n_theta=360;
    n_r=ceil(size(average,1)*0.5);
    [pf,~]=cryo_pft(masked_average,n_r,n_theta,'single');  % take Fourier transform of projections
    
    % Find common lines from projections
    log_message('Computeing common lines');
    max_shift=15;
    shift_step=1;
    [clstack,~,~,~]=...
        cryo_clmatrix_gpu(pf,K,1,max_shift,shift_step);
    
    log_message('Saving common lines');
    matname=sprintf('abinitio_info_nn%d_nm%d_group%d',nnavg,nmeans,groupid);
    save(fullfile(workflow.info.working_dir,matname),...
        'n_theta','n_r','clstack','max_shift','shift_step');
    
    log_message('Starting buliding synchronization matrix');
    S=cryo_syncmatrix_vote(clstack,n_theta);
    log_message('Finished buliding synchronization matrix');
    
    rotations=cryo_syncrotations(S);
    
    save(fullfile(workflow.info.working_dir,matname),...
        'rotations','S','-append');
    
    d=eigs(S,10);
    d_str=sprintf('%7.2f ',d);
    log_message('Top 10 eigenvalues of sync matrix are %s',d_str);
    
    clerr=cryo_syncconsistency(rotations,clstack,n_theta);
    h=figure;
    hist(clerr(:,3),360);
    clerrFIGname=fullfile(workflow.info.working_dir,...
        sprintf('clerr_nn%02d_group%d.fig',nnavg,groupid));
    clerrEPSname=fullfile(workflow.info.working_dir,...
        sprintf('clerr_nn%02d_group%d.eps',nnavg,groupid));
    hgsave(clerrFIGname);
    print('-depsc',clerrEPSname);
    close(h);

    log_message('Estimating shifts');
    [est_shifts,~]=cryo_estimate_shifts(pf,rotations,max_shift,shift_step,10000,[],0);
    log_message('Finished estimating shifts');
    
    save(fullfile(workflow.info.working_dir,matname),'est_shifts','-append');
    
    % Reconstruct downsampled volume with no CTF correction
    n=size(average,1);
    [ v1, ~, ~ ,~, ~, ~] = recon3d_firm( average,...
        rotations,-est_shifts, 1e-6, 100, zeros(n,n,n));
    ii1=norm(imag(v1(:)))/norm(v1(:));
    log_message('Relative norm of imaginary components = %e\n',ii1);
    v1=real(v1);
    volname=sprintf('vol_nn%d_nm%d_group%d.mrc',nnavg,nmeans,groupid);
    WriteMRC(v1,1,fullfile(workflow.info.working_dir,volname));
    log_message('Saved %s',volname);
    
    % Reconstruct downsampled volume with CTF correction
    fname=sprintf('ctfs_effective_nn%02d_group%d.mrc',nnavg,groupid);
    fullfilename=fullfile(workflow.info.working_dir,fname);
    
    if exist(fullfilename,'file')~=2
        log_message('%s does not exist. Skipping reconstruction with CTF correction',fullfilename);
    else
        log_message('Reconstructing with CTF correction.');
        ctfs=ReadMRC(fullfilename);
        [ v1, ~, ~,~, ~, ~] = recon3d_firm_ctf( average,...
            ctfs(:,:,1:size(average,3)),1:size(average,3),rotations,-est_shifts,...
            1.0e-6, 100, zeros(n,n,n));
        ii1=norm(imag(v1(:)))/norm(v1(:));
        log_message('Relative norm of imaginary components = %e\n',ii1);
        v1=real(v1);
        volname=sprintf('vol_ctf_corrected_nn%d_nm%d_group%d.mrc',nnavg,nmeans,groupid);
        WriteMRC(v1,1,fullfile(workflow.info.working_dir,volname));
        log_message('Saved %s',volname);
    end
    
end

log_message('Workflow file: %s\n',workflow_fname);
log_message('Use this file name when calling subsequent funtions.\n');
log_message('Call next cryo_workflow_abinitio_analyze(''%s'')\n',workflow_fname);
