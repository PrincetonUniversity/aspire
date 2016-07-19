function PFprojs=cryo_phaseflip(CTFdata,projs,prefix)
% CRYO_PHASEFLIP Phase flip projections
%
% PFprojs=cryo_phaseflip(CTFdata,projs)   
%   Apply phase flipping to all projections in the stack projs. CTFdata
%   contains one record (with CTF parameters) for each projection. It is
%   generated, for example, by reading a STAR file (see example below).
%   CTF records with no corresponding projetctions are ignored.
%
% PFprojs=cryo_phaseflip(CTFdata,N,prefix) 
%   Read the projections from disk according to the image name given in
%   each CTF record. Only the first N records are processed. prefix is
%   added to image names as only relative filenames are stored in CTFdata.
% 
% PFprojs=cryo_phaseflip(CTFdata,projs,pixA) 
%   Use pixA for pixel size (in Angstroms).

% Example:
%   CTFdata=readSTAR(fname);
%   FPprojs=cryo_phaseflip(CTFdata,100); % Phase flip the first 100 projections
%
% Yoel Shkolnisky, July 2014.
%
% Revisions:
% Y.S. Novmber 2014     Revise function's interface.
% Y.S. March 2015       Add option to provide pixel size.
% Y.S. August 2015      Look for pixA in the STAR file.

stackgiven=1;
if isscalar(projs) % No stack is given, just number of projections.
    Nprojs=projs;
    projs=-1;
    stackgiven=0;
else
    Nprojs=size(projs,3);
    if nargin==3
        pixA=prefix;
    end
end

lastStackProcessed='';
N=numel(CTFdata.data);
MRCstack=-1;
projsinit=0; % Has the stack PFprojs been initialized already.
       
Nprojs=min(N,Nprojs);
%fprintf('Phaseflip:')
printProgressBarHeader;
for k=1:Nprojs
    progressTic(k,Nprojs);
    if ~stackgiven % No stack is give
        % Get the identification string of the next image to process.
        % This is composed from the index of the image within an image stack,
        % followed by '@' and followed by the filename of the MRC stack.
        imageID=CTFdata.data{k}.rlnImageName;
        imparts=strsplit(imageID,'@');
        imageidx=str2double(imparts{1});
        stackname=imparts{2};
        
        % Read the image stack from the disk, if different from the current
        % one.
        if ~strcmp(stackname,lastStackProcessed)
            MRCstack=ReadMRC(fullfile(prefix,stackname));
            lastStackProcessed=stackname;
        end
        
        if imageidx>size(MRCstack,3)
            error('image %d in stack %s does not exist',imageidx, stackname);
        end
        
        % Get CTF parameters for the given image.
        im=MRCstack(:,:,imageidx);
    else
        stackname='N/A';
        imageidx=k;
        im=projs(:,:,k);        
    end
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
    
    [voltage,DefocusU,DefocusV,DefocusAngle,Cs,pixA,A]=...
        cryo_parse_Relion_CTF_struct(CTFdata.data{k});
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
