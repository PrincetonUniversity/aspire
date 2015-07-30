function cryo_workflow_preprocess(workflow_fname)

% Read workflow file
tree=xmltree(workflow_fname);
workflow=convert(tree);

% Check that the input has not changed since the last time is was
% preprocessed. % Compute hash of data file and compare to the stored hash.
log_message('Validating input...\n');
opt.Input='file';
opt.Format='hex';
opt.Method='MD5';
% hash=DataHash(workflow.info.rawdata,opt);
% if ~strcmpi(hash,workflow.info.rawdatahash)
%     message='Raw projection''s file has changed. Proceed? ';
%     do_proceed=multichoice_question(message,{'Y','N'},[ 1, 0],'Y');
%     if ~do_proceed
%        log_message('Aborting...\n');
%        return;
%     end
% end

%% Read preprocessing parameters
% Do phase flip?
message='Phaseflip projections? ';
do_phaseflip=multichoice_question(message,{'Y','N'},[ 1, 0],'Y');
if do_phaseflip==1
    ctfdata=fmtinput('Enter full path of CTF file: ','','%s');
end

% Do downsample?
[~,mrc_header]=ReadMRC(workflow.info.rawdata,1, 1);
% Input file has  mrc_header.nz images, each of size mrc_header.nx x
% mrc_header.ny
fprintf('%s :%dx%d (%d projections)\n',...
    workflow.info.rawdata,mrc_header.nx,mrc_header.ny,mrc_header.nz);
message='Number of projections to read? ';
nprojs=fmtinput(message,mrc_header.nz,'%d');

message='Downsample? ';
do_downsample=multichoice_question(message,{'Y','N'},[ 1, 0],'Y');
if do_downsample==1
    message='Downsample to size? ';
    projdim=fmtinput(message,mrc_header.nx,'%d');
end

% Do normalize background.
message='Normalize background of images to variance 1? ';
do_normalize=multichoice_question(message,{'Y','N'},[ 1, 0],'Y');

% Do prewhiten
message='Prewhiten? ';
do_prewhiten=multichoice_question(message,{'Y','N'},[ 1, 0],'Y');

% Do split into groups?
message='Split data into groups? ';
do_split=multichoice_question(message,{'Y','N'},[ 1, 0],'Y');
if do_split==1
    numgroups=fmtinput('Number of groups ',2,'%d');
else
    numgroups=1;
end

%% Update workflow struct
workflow.preprocess.phaseflip=do_phaseflip;
workflow.preprocess.ctfdata=ctfdata;
workflow.preprocess.nprojs=nprojs;
workflow.preprocess.do_downsample=do_downsample;
workflow.preprocess.projdim=projdim;
workflow.preprocess.do_normalize=do_normalize;
workflow.preprocess.do_prewhiten=do_prewhiten;
workflow.preprocess.split=do_split;
workflow.preprocess.numgroups=numgroups;

tree=struct2xml(workflow);
save(tree,workflow_fname); 


%% Execute preprocessing

open_log(fullfile(workflow.info.working_dir,workflow.info.logfile));

% Load data
log_message('Loading data\n');
projs=ReadMRC(workflow.info.rawdata,1,nprojs);

% Phaseflip
if do_phaseflip
    log_message('Reading CTF data\n');
    CTFdata=readSTAR(workflow.preprocess.ctfdata);
    log_message('Phasdflipping\n');
    PFprojs=cryo_phaseflip(CTFdata,projs);
else
    log_message('Skipping phaseflip\n');
    PFprojs=projs;
end
clear projs

% Downsample
if do_downsample
    log_message('Downsampling\n');
    PFDprojs=cryo_downsample(PFprojs,[projdim projdim],1); 
else
    log_message('Skipping downsampling\n');
    PFDprojs=PFprojs;
end
clear PFprojs

% Normalize images
if do_normalize
    log_message('Normalize background\n');
    n=size(PFDprojs,1);
    PFDprojs=cryo_normalize_background(PFDprojs,round(n/2)-10);
end

% Prewhiten
if do_prewhiten
    % Estimate noise PSD and prewhiten
    log_message('Estimating noise power spectrum');     
    n=size(PFDprojs,1);
    log_message('Each projection of size %d x %d',n,n);
    psd = Noise_Estimation(PFDprojs);
    log_message('Finished noise power spectrum estimation');    
    
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
    
    psd_white=Noise_Estimation(prewhitened_projs);
    
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
numgroups=workflow.preprocess.numgroups;
K2=floor(K/numgroups);

for groupid=1:numgroups
    fname=sprintf('phaseflipped_downsampled_prewhitened_group%d.mrc',groupid);
    fullfilename=fullfile(workflow.info.working_dir,fname);
    log_message('Saving group %d\n',groupid);
    WriteMRC(single(prewhitened_projs(:,:,(groupid-1)*K2+1:groupid*K2)),1,fullfilename);
end
clear prewhitened_projs

xmlname=sprintf('%s.xml',workflow_fname);
log_message('Workflow file: %s\n',xmlname);
log_message('Use this file name when calling subsequent funtions.\n');
log_message('Call next cryo_workflow_classify(''%s'')\n',xmlname);

close_log;