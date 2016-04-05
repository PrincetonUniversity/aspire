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

matname=fullfile(workflow.info.working_dir,'preprocess_info'); % mat file
    % to save intermediate data.

% Load data
nprojs=str2double(workflow.preprocess.nprojs);
log_message('Loading data %d projections from %s',nprojs,workflow.info.rawdata);
projs=ReadMRC(workflow.info.rawdata,1,nprojs);
szprojs=size(projs);
if szprojs(1)~=szprojs(2)
    error('Input projections must be square.')
end


if ~exist(workflow.preprocess.ctfdata,'file')
    error('Cannot read STAR file %s',workflow.preprocess.ctfdata);
end

% Phaseflip
if str2double(workflow.preprocess.phaseflip)
    log_message('Reading CTF data %s',workflow.preprocess.ctfdata);
    CTFdata=readSTAR(workflow.preprocess.ctfdata);
    log_message('Phaseflipping');
    PFprojs=cryo_phaseflip(CTFdata,projs);
else
    log_message('Skipping phaseflip');
    PFprojs=projs;
end
clear projs

% Crop
if str2double(workflow.preprocess.do_crop)
    croppeddim=str2double(workflow.preprocess.croppeddim);
    log_message('Cropping to %dx%d',croppeddim, croppeddim);
    PFCprojs=cryo_crop(PFprojs,[croppeddim croppeddim],1); 
else
    log_message('Skipping cropping');
    PFCprojs=PFprojs;
end
clear FPprojs

% Downsample
pixelscaling=1; % How much pixel size is changed due to downsampling.
if str2double(workflow.preprocess.do_downsample)
    downsampleddim=str2double(workflow.preprocess.downsampleddim);
    log_message('Downsampling to %dx%d',downsampleddim, downsampleddim);
    PFDprojs=cryo_downsample(PFCprojs,[downsampleddim downsampleddim],1); 
    pixelscaling=size(PFCprojs,1)/downsampleddim;
else
    log_message('Skipping downsampling');
    PFDprojs=PFCprojs;
end
clear PFCprojs

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
shuffleidx=1:K;
if fieldexist(workflow,'preprocess','do_shuffle') && ...
        str2double(workflow.preprocess.do_shuffle)==1
    log_message('Shuffling data');
    shuffleidx=randperm(K);
end

% Save the indcies of the images after shuffling.
if exist(matname,'file')
    save(matname,'shuffleidx','-append');
else
    save(matname,'shuffleidx'); 
end

numgroups=str2double(workflow.preprocess.numgroups);
K2=floor(K/numgroups);

for groupid=1:numgroups
    fname=sprintf('phaseflipped_cropped_downsampled_prewhitened_group%d.mrc',groupid);
    fullfilename=fullfile(workflow.info.working_dir,fname);
    log_message('Saving group %d',groupid);
    WriteMRC(single(prewhitened_projs(:,:,shuffleidx((groupid-1)*K2+1:groupid*K2))),1,fullfilename);
    
    % Write indices of the raw images in this group
    fname=sprintf('raw_indices_group%d.mrc',groupid);
    fullfilename=fullfile(workflow.info.working_dir,fname);
    fid=fopen(fullfilename,'w');
    if fid<0
        error('Failed to open %s',fullfilename);
    end
    fprintf(fid,'%d\n',shuffleidx((groupid-1)*K2+1:groupid*K2));
    fclose(fid);
    
    
    % Generate CTF data for the downsampled images in STAR format.
    % XXX This code is ratherslow due to the cell arrays - optimize.
    if ~isempty(workflow.preprocess.ctfdata)
        log_message('Generating CTF for downsampled images');
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
        log_message('Pixels size of original images %d Angstroms',pixA);
        log_message('Pixels size of downsampled images %d Angstroms',pixAdownsampled);
        
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
        writeSTAR(CTFdownsampled,fullfilename);
        log_message('Finished saving CTF data of downsampled images into %s',fullfilename);
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
        writeSTAR(CTFraw,fullfilename);
        
    end
end
clear prewhitened_projs

log_message('Workflow file: %s',workflow_fname);
log_message('Use this file name when calling subsequent funtions.');
log_message('Call next cryo_workflow_classify(''%s'')',workflow_fname);

close_log;
