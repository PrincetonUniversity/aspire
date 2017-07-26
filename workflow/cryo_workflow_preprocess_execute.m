function cryo_workflow_preprocess_execute(workflow_fname)
% CRYO_WORKFLOW_PREPROCESS_EXECUTE  Execute data set preprocessing
%
% cryo_workflow_preprocess_execute(workflow_fname)
%   Preprocess the projections data set according to the parameters stored
%   in the file workflow_fname.
%
% See also cryo_workflow_preprocess
%
% Yoel Shkolnisky, August 2015.

%% Validate workflow file
cryo_workflow_preprocess_validate(workflow_fname);

%% Read workflow file
tree=xmltree(workflow_fname);
workflow=convert(tree);

%% Execute preprocessing
open_log(fullfile(workflow.info.working_dir,workflow.info.logfile));

log_message('Starting cryo_workflow_preprocess_execute');
log_message('Loaded XML file %s (MD5: %s)',workflow_fname,MD5(workflow_fname));

matname=fullfile(workflow.info.working_dir,'preprocess_info.mat'); % mat file
    % to save intermediate data.

% Load data
%nprojs=str2double(workflow.preprocess.nprojs);
%log_message('Loading data %d projections from %s',nprojs,workflow.info.rawdata);
%projs=ReadMRC(workflow.info.rawdata,1,nprojs);
%szprojs=size(projs);

% Create backup of the raw data
log_message('Computing MD5 of %s',workflow.info.rawdata);
log_message('Using raw data %s (MD5: %s)',...
    workflow.info.rawdata,MD5(workflow.info.rawdata));
log_message('Creating backup of raw data.');
SRCname=tempmrcname; % SRC stands for "source".
log_message('Copying raw data from %s to temporary file %s',workflow.info.rawdata,SRCname);
copyfile(workflow.info.rawdata,SRCname);

rawdataReader=imagestackReader(SRCname);
szprojs=rawdataReader.dim;
if szprojs(1)~=szprojs(2)
    error('Input projections must be square.')
end

% Phaseflip
if str2double(workflow.preprocess.phaseflip)
    log_message('Start phaseflipping...')
    if ~exist(workflow.preprocess.ctfdata,'file')
        error('Cannot read STAR file %s',workflow.preprocess.ctfdata);
    end
    log_message('Reading CTF data %s (MD5: %s)',...
        workflow.preprocess.ctfdata,...
        MD5(workflow.preprocess.ctfdata));
    
    CTFdata=readSTAR(workflow.preprocess.ctfdata);
    %PFprojs=cryo_phaseflip(CTFdata,projs);
    PFfname=tempmrcname; %PF stands for phaseflipped
    log_message('Phaseflipped images will be saved to temporary file %s',PFfname);
    log_message('Running cryo_phaseflip_outofcore');
    cryo_phaseflip_outofcore(CTFdata,SRCname,PFfname);
    log_message('Finished phaseflipping')
    
    delete(SRCname);
else
    log_message('Skipping phaseflip');
    % Copy raw images to temporary file
    PFfname=SRCname; % No phaseflip so continute with the raw data file.
end

%PFprojs=ReadMRC(PFfname);
%clear projs

% Crop
if str2double(workflow.preprocess.do_crop)
    croppeddim=str2double(workflow.preprocess.croppeddim);
    log_message('Start cropping...');
    log_message('Cropping to %dx%d',croppeddim, croppeddim);
    % PFCprojs=cryo_crop(PFprojs,[croppeddim croppeddim],1); 
    PFCfname=tempmrcname; % PFC stands for phaseflipped+cropped
    log_message('Cropped images will be saved to temporary file %s',PFCfname);
    log_message('Running cryo_crop_outofcore');
    cryo_crop_outofcore(PFfname,PFCfname,[croppeddim croppeddim])
    log_message('Finished cropping');    
    
    delete(PFfname);
else
    log_message('Skipping cropping');
    %PFCprojs=PFprojs;
    PFCfname=PFfname;
end
%clear FPprojs

% Downsample
pixelscaling=1; % How much pixel size is changed due to downsampling.
if str2double(workflow.preprocess.do_downsample)
    log_message('Start downsampling...');
    downsampleddim=str2double(workflow.preprocess.downsampleddim);
    log_message('Downsampling to %dx%d',downsampleddim, downsampleddim);
    %PFDprojs=cryo_downsample(PFCprojs,[downsampleddim downsampleddim],1); 
    PFCDfname=tempmrcname; %PFCD stands for phaseflipped+cropped+downsampled
    log_message('Downsampled images will be saved to temporary file %s',PFCDfname);
    log_message('Running cryo_downsample_outofcore');
    cryo_downsample_outofcore(PFCfname,PFCDfname,[downsampleddim downsampleddim]);
    log_message('Finished downsampling');
    
    PFCReader=imagestackReader(PFCfname);
    originaldim=PFCReader.dim(1);
    pixelscaling=originaldim/downsampleddim;
    
    delete(PFCfname)
else
    log_message('Skipping downsampling');
    %PFDprojs=PFCprojs;
    PFCDfname=PFCfname;
end
%clear PFCprojs

% Normalize images
if str2double(workflow.preprocess.do_normalize)
    log_message('Start normalize background...');
    % n=size(PFDprojs,1);
    % PFDprojs=cryo_normalize_background(PFDprojs,round(n/2)-10);

    PFCDReader=imagestackReader(PFCDfname);
    n=PFCDReader.dim(1);
    PFCDNfname=tempmrcname; % phaseflipped+cropped+downsampled+normalized
    log_message('Normalized images will be saved to temporary file %s',PFCDNfname);
    log_message('Running cryo_normalize_background_outofcore');
    cryo_normalize_background_outofcore(PFCDfname,PFCDNfname,round(n/2)-10);
    log_message('Finished normalizing');
    
    delete(PFCDfname);
else
    log_message('Skipping background normaliztion');
    PFCDNfname=PFCDfname;
end

% Prewhiten
if str2double(workflow.preprocess.do_prewhiten)
    log_message('Starting prewhitening...');
    % Estimate noise PSD and prewhiten
    log_message('Estimating noise power spectrum');     
    % n=size(PFDprojs,1);
    PFCDNReader=imagestackReader(PFCDNfname);
    n=PFCDNReader.dim(1);
    log_message('Each projection of size %dx%d',n,n);
    %psd = cryo_noise_estimation(PFDprojs);
    psd = cryo_noise_estimation_outofcore(PFCDNfname);
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
    
    
    %log_message('Prewhitening images');
    %prewhitened_projs = cryo_prewhiten(PFDprojs, psd);
    %fname=sprintf('phaseflipped_downsampled_prewhitened.mrc');
    %WriteMRC(single(prewhitened_projs),1,fullfile(workflow.info.working_dir,fname));

    PFCDNWfname=tempmrcname;
        % phaseflipped+cropped+downsampled+normalized+whitened
    log_message('Prewhitened images will be saved to temporary file %s',PFCDNWfname);
    log_message('Running cryo_prewhiten_outofcore');
    cryo_prewhiten_outofcore(PFCDNfname,PFCDNWfname,psd);
    log_message('Finished prewhitening');
    
    % Compute power spectrum of the prewhitened images - just to verification
    % Normalize projections to norm 1
    log_message('Computing power spectrum of prewhitened projections - for verifying that power spectrum is white');
    
    %psd_white=cryo_noise_estimation(prewhitened_projs);
    psd_white = cryo_noise_estimation_outofcore(PFCDNWfname);
    
    h=figure;
    plot(psd_white(n,:));
    title('Noise spectrum of prewhitened-projections');
    psdFIGname=fullfile(workflow.info.working_dir,'psd_after_prewhitening.fig');
    psdEPSname=fullfile(workflow.info.working_dir,'psd_after_prewhitening.eps');
    hgsave(psdFIGname);
    print('-depsc',psdEPSname);
    close(h);
    
    delete(PFCDNfname);

else
    PFCDNWfname=PFCDNfname;
end

%clear PFDprojs

% Global phase flip
log_message('Starting global phaseflip...') 
%[prewhitened_projs,doflip]=cryo_globalphaseflip(prewhitened_projs);
fname=sprintf('phaseflipped_downsampled_prewhitened.mrc');
PFCDNWGfname=fullfile(workflow.info.working_dir,fname);
doflip=cryo_globalphaseflip_outofcore(PFCDNWfname ,PFCDNWGfname);

if doflip
    log_message('Phase of images was flipped');
else
    log_message('No need to global phase flip. Phase of images not flipped');
end
log_message('Images copied (even if not globally phase flipped) to %s (MD5: %s)'...
    ,PFCDNWGfname,MD5(PFCDNWGfname));
log_message('Finished global phaseflip...')

delete(PFCDNWfname);

% Split into groups
PFCDNWGReader=imagestackReader(PFCDNWGfname,1); % Cache size is 1 since we 
        % expect to shuffle the images which would result in many cache
        % flush.
K=PFCDNWGReader.dim(3);
%K=size(prewhitened_projs,3);
shuffleidx=1:K;
if fieldexist(workflow,'preprocess','do_shuffle') && ...
        str2double(workflow.preprocess.do_shuffle)==1
    log_message('Images will be shuffled prior to splitting to groups');
    log_message('Shuffled order of images saved to variable ''shuffleidx'' in %s',matname);
    shuffleidx=randperm(K);
else
    log_message('Images will NOT be shuffled prior to splitting to groups');
end

% Save the indcies of the images after shuffling.
if exist(matname,'file')
    save(matname,'shuffleidx','-append');
else
    save(matname,'shuffleidx'); 
end

numgroups=str2double(workflow.preprocess.numgroups);
K2=floor(K/numgroups);

log_message('Start splitting to %d groups',numgroups);
for groupid=1:numgroups
    fname=sprintf('phaseflipped_cropped_downsampled_prewhitened_group%d.mrc',groupid);
    fullfilename=fullfile(workflow.info.working_dir,fname);
    log_message('Saving group %d into file %s',groupid,fullfilename);
    %WriteMRC(single(prewhitened_projs(:,:,shuffleidx((groupid-1)*K2+1:groupid*K2))),1,fullfilename);
    
    groupstack=imagestackWriter(fullfilename,K2,1);
    for k=1:K2
        proj=PFCDNWGReader.getImage(shuffleidx((groupid-1)*K2+k));
        groupstack.append(proj);
    end
    groupstack.close;    
    log_message('Saved group %d into file %s (MD5: %s)',...
        groupid,fullfilename,MD5(fullfilename));
    
    % Write indices of the raw images in this group
    fname=sprintf('raw_indices_group%d.dat',groupid);
    fullfilename=fullfile(workflow.info.working_dir,fname);    
    fid=fopen(fullfilename,'w');
    if fid<0
        error('Failed to open %s',fullfilename);
    end
    fprintf(fid,'%d\n',shuffleidx((groupid-1)*K2+1:groupid*K2));
    fclose(fid);
    log_message('Indices of raw images in group %d are saved to %s (MD5: %s)'...
        ,groupid,fullfilename,MD5(fullfilename));
    
    
    % Generate CTF data for the downsampled images in STAR format.
    % XXX This code is ratherslow due to the cell arrays - optimize.
    if ~isempty(workflow.preprocess.ctfdata)
        log_message('Generating CTF for downsampled images of group %d',groupid);
        CTFdownsampled=createSTARdata(K2,'rlnVoltage','rlnDefocusU','rlnDefocusV',...
            'rlnDefocusAngle','rlnSphericalAberration',...
            'rlnAmplitudeContrast','pixA','phaseflipped');
        
        
        % Compute pixel size of downsampled images.
        idx=shuffleidx((groupid-1)*K2+1); % Read parameters of the first image
        % in the current group.
        [~,~,~,~,~,pixA,~]=cryo_parse_Relion_CTF_struct(CTFdata.data{idx});
        % Read pixel size
        pixAdownsampled=pixA*pixelscaling;
        workflow.preprocess.pixA=pixA;
        workflow.preprocess.pixAdownsampled=pixAdownsampled;
        tree=struct2xml(workflow);
        save(tree,workflow_fname);
        log_message('Pixel size of original images %4.2f Angstroms',pixA);
        log_message('Pixel size of downsampled images %4.2f Angstroms',pixAdownsampled);
        
        log_message('Create CTF data for downsampled images');
        printProgressBarHeader;
        
        for k=1:K2
            progressTicFor(k,K2);
            idx=shuffleidx((groupid-1)*K2+k);
            [voltage,DefocusU,DefocusV,DefocusAngle,Cs,pixA,A]=...
                cryo_parse_Relion_CTF_struct(CTFdata.data{idx});
            % Convert defocus to angstroms and defocus angle to degrees.
            CTFdownsampled=addrecordtoSTARdata(CTFdownsampled,k,'rlnVoltage',voltage,...
                'rlnDefocusU',DefocusU*10,'rlnDefocusV',DefocusV*10,...
                'rlnDefocusAngle',DefocusAngle*180/pi,...
                'rlnSphericalAberration',Cs,...
                'rlnAmplitudeContrast',A,'pixA',pixA*pixelscaling,...
                'phaseflipped',str2double(workflow.preprocess.phaseflip));
            
            if workflow.preprocess.pixA~=pixA
                error('Pixel size varies in original data');
            end
        end
        
        fname=sprintf('ctfs_group%d.star',groupid);
        fullfilename=fullfile(workflow.info.working_dir,fname);
        log_message('Saving CTF data for group %d to %s',groupid,fullfilename);
        writeSTAR(CTFdownsampled,fullfilename);
        %log_message('Finished saving CTF data of downsampled images into %s',fullfilename);

        %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %     % Generate filter reponse of the CTFs of the downsampled images.
        %     % Can be used for testing the above code.
        %
        %     log_message('Generating CTF for group %d',groupid);
        %     n=size(prewhitened_projs,1);
        %     CTFds=zeros(n,n,K2,class(prewhitened_projs)); % CTF response
        %     % of the downsampled images in each group.
        %     for k=1:K2
        %         idx=shuffleidx((groupid-1)*K2+k);
        %         [voltage,DefocusU,DefocusV,DefocusAngle,Cs,pixA,A]=...
        %             cryo_parse_Relion_CTF_struct(CTFdata.data{idx});
        %         h=cryo_CTF_Relion(n,voltage,DefocusU,DefocusV,DefocusAngle,...
        %             Cs,pixA*pixelscaling,A);
        %         if str2double(workflow.preprocess.phaseflip)
        %             h=abs(h);
        %         end
        %         CTFds(:,:,k)=h;
        %     end
        %
        %     fname=sprintf('ctfs_group%d.mrc',groupid);
        %     fullfilename=fullfile(workflow.info.working_dir,fname);
        %     WriteMRC(single(CTFds),1,fullfilename);
        %
        %     CTFdownsampled_response=zeros(n,n,K2,class(prewhitened_projs)); % CTF response
        %
        %     for k=1:K2
        %         [voltage,DefocusU,DefocusV,DefocusAngle,Cs,pixA,A]=...
        %             cryo_parse_Relion_CTF_struct(CTFdownsampled.data{k});
        %         h=cryo_CTF_Relion(n,voltage,DefocusU,DefocusV,DefocusAngle,...
        %             Cs,pixA,A);
        %         if str2double(workflow.preprocess.phaseflip)
        %             h=abs(h);
        %         end
        %         CTFdownsampled_response(:,:,k)=h;
        %     end
        %
        %     assert(norm(CTFdownsampled_response(:)-CTFds(:))/norm(CTFds(:))<1.0e-14);
        %
        %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Save CTF parameters of raw images in each group. We don't save the
        % reponse for each image since these images may be very large.
        
        CTFraw=CTFdata;
        CTFraw.data=CTFraw.data(shuffleidx((groupid-1)*K2+1:groupid*K2));
        fname=sprintf('ctfs_raw_group%d.star',groupid);
        fullfilename=fullfile(workflow.info.working_dir,fname);
        log_message('Save CTF parameters of raw images group %d into %s',groupid,fullfilename);
        writeSTAR(CTFraw,fullfilename);
        
    end
end

delete(PFCDNWGfname);

log_message('Finished splitting to groups');
%clear prewhitened_projs

log_message('Data save into MAT file %s (MD5: %s)',matname,MD5(matname));

log_message('Workflow file: %s',workflow_fname);
log_message('Use this file name when calling subsequent funtions.');
log_message('Call next cryo_workflow_classify(''%s'')',workflow_fname);

close_log;
