function cryo_phaseflip_outofcore(CTFdata,prefix,outstackname,pixA,verbose)
% CRYO_PHASEFLIP_OUTOFCORE      Phase flip projections
%
% cryo_phaseflip_outofcore(CTFdata,prefix,outstackname)
%   Apply phase flipping to all projections specified in CTFdata, and write
%   the flipped images to the MRC file outstackname.
%   CTFdata  contains one record (with CTF parameters) for each projection
%   in the input stack. It is generated, for example, by reading a STAR
%   file (see example below). 
%   prefix is added to image names as only relative filenames are stored in
%   CTFdata. Pass prefix as empty string of [] if prefix is not needed.
%   The function does not load the entire stack into memory.
%   If pixA==-1, then the function is trying to read pixel size from the
%   CTFdata. Otherwise, uses the give pixA.
%
% Example:
%   CTFdata=readSTAR(fname);
%   cryo_phaseflip_outofcore(CTFdata,'outstack.mrcs'); 
%
% Yoel Shkolnisky, May 2016.
%
% Revisions:
% Y.S. October 2017     Eliminate input stack name. Images are always read
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
% NOTE: If STAR files of RELION 3.1 is used, then the structure of the
% STAR file is assumed to contained one optics group (location 1 in the
% stardata array) and one particles group (location 2 in the stardata
% array).
if numel(CTFdata)==1 % RELION version < 3.1
    Nprojs=numel(CTFdata.data);
else
    Nprojs=numel(CTFdata(2).data);
end
outstack=imagestackWriter(outstackname,Nprojs,pixA);

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
    if numel(CTFdata)==1 % RELION version < 3.1
        imageID=CTFdata.data{k}.rlnImageName;
    else
        imageID=CTFdata(2).data{k}.rlnImageName;
    end
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
        MRCstack=imagestackReader(MRCname);
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

    [voltage,DefocusU,DefocusV,DefocusAngle,Cs,tmppixA,A]=...
        cryo_parse_Relion_CTF_struct(CTFdata,k);
    
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
        % This test is only valid for odd n.
        if norm(imag(pfim(:)))/norm(pfim(:))>5.0e-7 % images are single precision.
            warning('Large imaginary components in image %d = %e',imageidx,norm(imag(pfim(:)))/norm(pfim(:)));
        end
    end
    pfim=real(pfim);
    outstack.append(single(pfim));
    
    if verbose==2
        if mod(k,1000)==0
            log_message('Finsihed processing %d/%d projections',k,Nprojs);
        end
    end    
end
outstack.close;
