function PFprojs=cryo_phaseflip(CTFdata,prefix,pixA,verbose)
% CRYO_PHASEFLIP Phase flip projections
%
% PFprojs=cryo_phaseflip(CTFdata,prefix) 
%   Read the projections from disk according to the image name given in
%   each CTF record. prefix is added to image names as only relative
%   filenames are stored in CTFdata.
% 
% PFprojs=cryo_phaseflip(CTFdata,prefix,pixA) 
%   Use pixA for pixel size (in Angstroms).
%   If pixA==-1, then the function is trying to read pixel size from the
%   CTFdata. Otherwise, uses the give pixA.
%   Pass prefix as empty string of [] if prefix is not needed.
%
% PFprojs=cryo_phaseflip(CTFdata,prefix,pixA,verbose) 
%   Set verbose to 1 to print progress bar 
%   set verbose to 2 to summarized progress info.
%   Set verbose to 3 to print messages on each processed projection.
%
%
% Example:
%   CTFdata=readSTAR(fname);
%   FPprojs=cryo_phaseflip(CTFdata); 
%
% Yoel Shkolnisky, July 2014.
%
% Revisions:
% Y.S. Novmber 2014     Revise function's interface.
% Y.S. March 2015       Add option to provide pixel size.
% Y.S. August 2015      Look for pixA in the STAR file.
% Y.S. October 2017     Eliminate input stack. Images are always read
%                       according to the name given in CTFdata

if ~exist('prefix','var')
    prefix='';
end

if isempty(prefix)
    prefix='';
end

if ~exist('verbose','var')
    verbose=0;
end

lastStackProcessed='';
Nprojs=numel(CTFdata.data);
projsinit=0; % Has the stack PFprojs been initialized already.

if verbose==1
    printProgressBarHeader;
end
for k=1:Nprojs
    if verbose==1
        progressTic(k,Nprojs);
    end
    % Get the identification string of the next image to process.
    % This is composed from the index of the image within an image stack,
    % followed by '@' and followed by the filename of the MRC stack.
    imageID=CTFdata.data{k}.rlnImageName;
    imparts=strsplit(imageID,'@');
    imageidx=str2double(imparts{1});
    stackname=imparts{2};
    
    if verbose>=3
        log_message('Processing prjection %d/%d: reading image %d from stack %s',...
            k,Nprojs,imageidx,stackname);
    end
    
    % Read the image stack from the disk, if different from the current
    % one.
    if ~strcmp(stackname,lastStackProcessed)
        MRCname=fullfile(prefix,stackname);
        MRCstack=imagestackReader(MRCname,100);
        if verbose>=2
            log_message('Processing prjection %d/%d: Initializing MRC stack %s',...
                k,Nprojs,MRCname);
        end
        lastStackProcessed=stackname;
    end
    
    if imageidx>MRCstack.dim(3)
        error('image %d in stack %s does not exist',imageidx, stackname);
    end
    
    % Get CTF parameters for the given image.
    im=MRCstack.getImage(imageidx);
    im=double(im); % Convert to double to eliminate small numerical
        % roundoff errors when comparing the current function to
        % cryo_phaseflip_outofcore. This line is only required to get
        % perfectly zero error when comparing the functions.
    n=size(im,1);
    if n~=size(im,2)
        error('Images must be square');
    end
    
    if projsinit==0
        PFprojs=zeros(n,n,Nprojs,'single');
        projsinit=1;
    end;
    
    [voltage,DefocusU,DefocusV,DefocusAngle,Cs,tmppixA,A]=...
        cryo_parse_Relion_CTF_struct(CTFdata.data{k});
    
    if verbose>=3
        log_message('Processing prjection %d/%d: CTF params read voltage=%d, DefocusU=%d, DefocusV=%d, DefocusAngle=%d, Cs=%d, pixA=%d ,A=%d',...
            k,Nprojs,voltage,DefocusU,DefocusV,DefocusAngle,Cs,tmppixA,A);
    end
    if pixA==-1
        if tmppixA~=-1
            pixA=tmppixA; % Use pixel size from STAR file
        else
            errmsg=sprintf(strcat('Pixel size not provided and does not appear in STAR file.\n',...
            'Provide it manually, or in the STAR file using the fields ',...
            'pixA or rlnDetectorPixelSize together with rlnMagnification.\n',...
            'You can use addfieldtoSTARdata to add pixA manually to the STAR file.'));
        error('%s',errmsg);
        end
    end
    
    h=cryo_CTF_Relion(n,voltage,DefocusU,DefocusV,DefocusAngle,Cs,pixA,A);
    
    % Phase flip
    imhat=fftshift(fft2(im));
    pfim=ifft2(ifftshift(imhat.*sign(h)));
    
    if mod(n,2)==1
        % This test is only vali for odd n.
        if norm(imag(pfim(:)))/norm(pfim(:))>5.0e-7 % images are single precision.
            warning('Large imaginary components in image %d in stack %s = %e',imageidx, stackname,norm(imag(pfim(:)))/norm(pfim(:)));
        end
    end
    pfim=real(pfim);
    PFprojs(:,:,k)=single(pfim);
    
    if verbose==2
        if mod(k,1000)==0
            log_message('Finsihed processing %d/%d projections',k,Nprojs);
        end
    end
    
%     clf;
%     subplot(1,3,1);
%     imagesc(im); axis off; axis image; colormap(gray); title('Original');
%     subplot(1,3,2);
%     imagesc(h); axis off; axis image; colormap(gray); title('CTF');
%     subplot(1,3,3);
%     imagesc(pfim); axis off; axis image; colormap(gray); title('Phase flipped');
%     fprintf('Processing image %d/%d from %s\n',imageidx,size(MRCstack,3),stackname);
%     pause(0.5);
end
