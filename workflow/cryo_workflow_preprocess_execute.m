function cryo_workflow_preprocess_execute(workflow_fname)

%% Validate workflow file
cryo_workflow_preprocess_validate(workflow_fname);

%% Read workflow file
tree=xmltree(workflow_fname);
workflow=convert(tree);

%% Execute preprocessing
open_log(fullfile(workflow.info.working_dir,workflow.info.logfile));

matname=fullfile(workflow.info.working_dir,'preprocess_info'); % mat file
    % to save intermediate data.

% Load data
nprojs=str2double(workflow.preprocess.nprojs);
log_message('Loading data %d projections from %s',nprojs,workflow.info.rawdata);
projs=ReadMRC(workflow.info.rawdata,1,nprojs);

% Phaseflip
if str2double(workflow.preprocess.phaseflip)
    log_message('Reading CTF data %s',workflow.preprocess.ctfdata);
    CTFdata=readSTAR(workflow.preprocess.ctfdata);
    log_message('Phasdflipping');
    PFprojs=cryo_phaseflip(CTFdata,projs);
else
    log_message('Skipping phaseflip');
    PFprojs=projs;
end
clear projs

% Downsample
if str2double(workflow.preprocess.do_downsample)
    projdim=str2double(workflow.preprocess.projdim);
    log_message('Downsampling to %dx%d',projdim, projdim);
    PFDprojs=cryo_downsample(PFprojs,[projdim projdim],1); 
else
    log_message('Skipping downsampling');
    PFDprojs=PFprojs;
end
clear PFprojs

% Normalize images
if str2double(workflow.preprocess.do_normalize)
    log_message('Normalize background');
    n=size(PFDprojs,1);
    PFDprojs=cryo_normalize_background(PFDprojs,round(n/2)-10);
end

% Prewhiten
if str2double(workflow.preprocess.do_prewhiten)
    % Estimate noise PSD and prewhiten
    log_message('Estimating noise power spectrum');     
    n=size(PFDprojs,1);
    log_message('Each projection of size %dx%d',n,n);
    psd = cryo_noise_estimation(PFDprojs);
    log_message('Finished noise power spectrum estimation');    

    save(matname,'psd'); % Save the estimated PSD.
    
    h=figure;
    plot(psd(n,:));
    title('Noise spectrum of raw projections');
    psdFIGname=fullfile(workflow.info.working_dir,'psd_before_prewhitening.fig');
    psdEPSname=fullfile(workflow.info.working_dir,'psd_before_prewhitening.eps');
    hgsave(psdFIGname);
    print('-depsc',psdEPSname);
    close(h);
    
    
    log_message('Prewhitening images');
    prewhitened_projs = Prewhiten_image2d(PFDprojs, psd);
    fname=sprintf('phaseflipped_downsampled_prewhitened.mrc');
    WriteMRC(single(prewhitened_projs),1,fullfile(workflow.info.working_dir,fname));
    log_message('Finished prewhitening images');
    
    % Compute power spectrum of the prewhitened images - just to verification
    % Normalize projections to norm 1
    log_message('Compute power spectrum of prewhitened projections - for verifying that power spectrum is white');
    
    psd_white=cryo_noise_estimation(prewhitened_projs);
    
    h=figure;
    plot(psd_white(n,:));
    title('Noise spectrum of prewhitened-projections');
    psdFIGname=fullfile(workflow.info.working_dir,'psd_after_prewhitening.fig');
    psdEPSname=fullfile(workflow.info.working_dir,'psd_after_prewhitening.eps');
    hgsave(psdFIGname);
    print('-depsc',psdEPSname);
    close(h);

else
    prewhitened_projs=PFDprojs;
end

clear PFDprojs

% Global phase flip
[prewhitened_projs,doflip]=cryo_globalphaseflip(prewhitened_projs);
if doflip
    log_message('Applying global phase flip');
end


% Split into groups
K=size(prewhitened_projs,3);
numgroups=str2double(workflow.preprocess.numgroups);
K2=floor(K/numgroups);

for groupid=1:numgroups
    fname=sprintf('phaseflipped_downsampled_prewhitened_group%d.mrc',groupid);
    fullfilename=fullfile(workflow.info.working_dir,fname);
    log_message('Saving group %d',groupid);
    WriteMRC(single(prewhitened_projs(:,:,(groupid-1)*K2+1:groupid*K2)),1,fullfilename);
end
clear prewhitened_projs

log_message('Workflow file: %s',workflow_fname);
log_message('Use this file name when calling subsequent funtions.');
log_message('Call next cryo_workflow_classify(''%s'')',workflow_fname);

close_log;
