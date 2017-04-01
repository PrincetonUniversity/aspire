function cryo_workflow_abinitio_execute(workflow_fname)
% CRYO_WORKFLOW_ABINITO_EXECUTE  Reconstruct abinitio models
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

log_message('Starting cryo_workflow_abinitio_execute');
log_message('Loaded XML file %s (MD5: %s)',workflow_fname,MD5(workflow_fname));

numgroups=str2double(workflow.preprocess.numgroups);
nmeans=str2double(workflow.abinitio.nmeans);
nnavg=str2double(workflow.abinitio.nnavg);

if isfield(workflow.abinitio,'algo')
    algo=str2double(workflow.abinitio.algo);
else
    algo=1; % Use sync3N by default
    log_message('workflow.abinitio.algo not found. Using sync3N for abinitio reconstruction.');
end

for groupid=1:numgroups
    reloadname=sprintf('averages_info_nn%02d_group%d.mat',nnavg,groupid);
    reloadname=fullfile(workflow.info.working_dir,reloadname);
    log_message('Loading %s (MD5: %s)',reloadname,MD5(reloadname));
    load(reloadname);
    fname=sprintf('averages_nn%02d_group%d.mrc',nnavg,groupid);
    fname=fullfile(workflow.info.working_dir,fname);
    log_message('Reading averages from %s (MD5: %s)',fname, MD5(fname));
    
    average=ReadMRC(fname);
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
        cryo_clmatrix(pf,K,1,max_shift,shift_step);
    
    log_message('Saving common lines');
    matname=sprintf('abinitio_info_nn%d_nm%d_group%d.mat',nnavg,nmeans,groupid);
    matname=fullfile(workflow.info.working_dir,matname);
    save(matname,'n_theta','n_r','clstack','max_shift','shift_step');
    
    % Orientation assigment
    if algo==1
        % Use sync3N
        VOTING_TICS_WIDTH=1;
        J_EIGS=4;
        J_WEIGHTS=true;
        S_WEIGHTS=true;
        
        % Estimate relative rotations
        [Rij0, r_valid_ks, r_good_ks, peak_width] = cryo_sync3n_estimate_all_Rijs...
            (clstack, n_theta, VOTING_TICS_WIDTH);
        
        log_message('Stage Done: Relative Rotations');
        
        % J-synchronization
        verbose = 2;
        [J_sync,J_significance,J_eigs,J_sync_iterations,~] = ...
            cryo_sync3n_Jsync_power_method (Rij0, J_EIGS, J_WEIGHTS, verbose);
        Rij = cryo_sync3n_flip_handedness(J_sync, Rij0);
        
        % Build 3NX3N Matrix S
        S = cryo_sync3n_syncmatrix(Rij);
        
        log_message('Stage Done: J Synchronization');
        
        % S Weighting
        if S_WEIGHTS
            % Estimate weights for the 3x3 blocks of S
            [W, Pij, scores_hist] = cryo_sync3n_syncmatrix_weights(Rij0);
        else
            W = ones(N); % Default weights are all one (no weights)
            Pij = [];
            scores_hist = struct();
        end
        
        log_message('Stage Done: Probabilistic Weighting');
        
        % Estimate rotations from S
        [rotations, S_eigs, ~] = cryo_sync3n_S_to_rot (S,10,W);
                
        d_str=sprintf('%7.2f ',S_eigs);
        log_message('Top 10 eigenvalues of (weighted) sync matrix are %s',d_str);
        
    elseif algo==2
        % Use sync2N
        log_message('Starting buliding synchronization matrix');
        S=cryo_syncmatrix_vote(clstack,n_theta);
        log_message('Finished buliding synchronization matrix');
        
        rotations=cryo_syncrotations(S);
                
        d=eigs(S,10);
        d_str=sprintf('%7.2f ',d);
        log_message('Top 10 eigenvalues of sync matrix are %s',d_str);
    elseif algo==3
        % Use LUD
        rotations = est_orientations_LUD(clstack, n_theta);
    else
        error('Invalid vaule for ''algo''.');
    end
    
    save(matname,'rotations','S','-append');

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
    
    save(matname,'est_shifts','-append');
    log_message('Saved %s (MD5: %s)',matname,MD5(matname));

    % Reconstruct downsampled volume with no CTF correction
    n=size(average,1);
    [ v1, ~, ~ ,~, ~, ~] = recon3d_firm( average,...
        rotations,-est_shifts, 1e-6, 100, zeros(n,n,n));
    ii1=norm(imag(v1(:)))/norm(v1(:));
    log_message('Relative norm of imaginary components = %e\n',ii1);
    v1=real(v1);
    volname=sprintf('vol_nn%d_nm%d_group%d.mrc',nnavg,nmeans,groupid);
    volname=fullfile(workflow.info.working_dir,volname);
    WriteMRC(v1,1,volname);
    log_message('Saved %s (MD5: %s)',volname,MD5(volname));
    
    % % % % %     % Reconstruct from the raw projections corresponding to each average.
    % The following code uses one raw projection for each class averages.
    % So if the number of averages is 1000, then only 1000 raw
    % images will be used.
    % % % % %     fname=sprintf('phaseflipped_cropped_downsampled_prewhitened_group%d.mrc',groupid);
    % % % % %     fullfilename=fullfile(workflow.info.working_dir,fname);
    % % % % %     log_message('Loading prewhitened projection from %s',fname);
    % % % % %     prewhitened_projs=ReadMRC(fullfilename);
    % % % % %
    % % % % %     reloadname=sprintf('averages_info_nn%02d_group%d',nnavg,groupid);
    % % % % %     reloadname=fullfile(workflow.info.working_dir,reloadname);
    % % % % %     alignment_data=load(reloadname);
    % % % % %
    % % % % %     % Note that at this point average is sorted according to increasing
    % % % % %     % norm_variance while prewithened_projs as well as all other variables
    % % % % %     % in the next call (except for average) are not.
    % % % % %     [ aligned_images ] = cryo_align_projections( prewhitened_projs,...
    % % % % %         alignment_data.VDM_angles, alignment_data.class_VDM,...
    % % % % %         alignment_data.class_VDM_refl, alignment_data.shifts,...
    % % % % %         nnavg, alignment_data.classcoreidx(1:1000),average,...
    % % % % %         alignment_data.doflip);
    % % % % %
    % % % % %
    % % % % %
    % % % % %     aligned_images2=zeros(size(aligned_images,1),size(aligned_images,2),...
    % % % % %         size(aligned_images,3)*size(aligned_images,4),'single');
    % % % % %     rotations2=zeros(3,3,size(aligned_images,3)*size(aligned_images,4));
    % % % % %     est_shifts2=zeros(size(aligned_images,3)*size(aligned_images,4),2,'single');
    % % % % %
    % % % % %     idx1=1;
    % % % % %     idx2=1;
    % % % % %     for k=1:size(aligned_images,3)
    % % % % %         for j=1:size(aligned_images,4)
    % % % % %             aligned_images2(:,:,idx1)=aligned_images(:,:,k,j);
    % % % % %             rotations2(:,:,idx1)=rotations(:,:,idx2);
    % % % % %             est_shifts2(idx1,:)=est_shifts(idx2,:);
    % % % % %             idx1=idx1+1;
    % % % % %         end
    % % % % %         idx2=idx2+1;
    % % % % %     end
    % % % % %
    % % % % %     [ v1, ~, ~ ,~, ~, ~] = recon3d_firm( aligned_images2,...
    % % % % %         rotations2,-est_shifts2, 1e-6, 100, zeros(n,n,n));
    % % % % %     ii1=norm(imag(v1(:)))/norm(v1(:));
    % % % % %     log_message('Relative norm of imaginary components = %e\n',ii1);
    % % % % %     v1=real(v1);
    % % % % %
    
%     % Reconstruct from the raw images used to generate the averages.
%     % All raw images of all class averages are used.
%     reloadname=sprintf('averages_info_nn%02d_group%d',nnavg,groupid);
%     reloadname=fullfile(workflow.info.working_dir,reloadname);
%     alignment_data=load(reloadname);
%     
%     averaging_data.rotations=rotations;
%     averaging_data.est_shifts=est_shifts;
%     averaging_data.nnavg=nnavg;
%     
%     log_message('Estimating shifts for raw projecitons');
%     [estR_raw,est_shifts_raw,rawimageindices]=cryo_assign_orientations_to_raw_projections(...
%         alignment_data,averaging_data);
%     
%     fname=sprintf('phaseflipped_cropped_downsampled_prewhitened_group%d.mrc',groupid);
%     fullfilename=fullfile(workflow.info.working_dir,fname);
%     log_message('Loading prewhitened projection from %s',fname);
%     prewhitened_projs=ReadMRC(fullfilename);
%     
%     n=size(average,1);
%     [ v1, ~, ~ ,~, ~, ~] = recon3d_firm(...
%         prewhitened_projs(:,:,rawimageindices),estR_raw,-est_shifts_raw,...
%         1e-6,100, zeros(n,n,n));
%     ii1=norm(imag(v1(:)))/norm(v1(:));
%     log_message('Relative norm of imaginary components = %e\n',ii1);
%     v1=real(v1);
%     volname=sprintf('vol_raw_nn%d_nm%d_group%d.mrc',nnavg,nmeans,groupid);
%     WriteMRC(v1,1,fullfile(workflow.info.working_dir,volname));
%     log_message('Saved %s',volname);
%     
%     save(fullfile(workflow.info.working_dir,matname),...
%         'estR_raw','est_shifts_raw','rawimageindices','-append');
    
    %     % Reconstruct downsampled volume with CTF correction
    %     fname=sprintf('ctfs_effective_nn%02d_group%d.mrc',nnavg,groupid);
    %     fullfilename=fullfile(workflow.info.working_dir,fname);
    %
    %     if exist(fullfilename,'file')~=2
    %         log_message('%s does not exist. Skipping reconstruction with CTF correction',fullfilename);
    %     else
    %         log_message('Reconstructing with CTF correction.');
    %         ctfs=ReadMRC(fullfilename);
    %         [ v1, ~, ~,~, ~, ~] = recon3d_firm_ctf( average,...
    %             ctfs(:,:,1:size(average,3)),1:size(average,3),rotations,-est_shifts,...
    %             1.0e-6, 100, zeros(n,n,n));
    %         ii1=norm(imag(v1(:)))/norm(v1(:));
    %         log_message('Relative norm of imaginary components = %e\n',ii1);
    %         v1=real(v1);
    %         volname=sprintf('vol_ctf_corrected_nn%d_nm%d_group%d.mrc',nnavg,nmeans,groupid);
    %         WriteMRC(v1,1,fullfile(workflow.info.working_dir,volname));
    %         log_message('Saved %s',volname);
    %      end
    
end

log_message('Workflow file: %s\n',workflow_fname);
log_message('Use this file name when calling subsequent funtions.\n');
log_message('Call next cryo_workflow_abinitio_analyze(''%s'')\n',workflow_fname);

close_log;