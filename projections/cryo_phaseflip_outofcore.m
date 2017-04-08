function cryo_phaseflip_outofcore(CTFdata,instackname,outstackname)
% CRYO_PHASEFLIP_OUTOFCORE      Phase flip projections
%
% cryo_phaseflip_outofcore(CTFdata,instackname,outstackname)
%   Apply phase flipping to all projections in the MRC file named
%   instackname, and write the flipped images to the MRC file outstackname.
%   CTFdata  contains one record (with CTF parameters) for each projection
%   in the input stack. It is generated, for example, by reading a STAR
%   file (see example below). 
%   The function does not load the entire stack into memory.
%
% Example:
%   CTFdata=readSTAR(fname);
%   cryo_phaseflip_outofcore(CTFdara,'instack.mrc','outstack.mrc'); 
%
% Yoel Shkolnisky, May 2016.


N=numel(CTFdata.data);
instack=imagestackReader(instackname,100);
Nprojs=min(N,instack.dim(3));
outstack=imagestackWriter(outstackname,Nprojs,1,100);

printProgressBarHeader;
for k=1:Nprojs
    progressTic(k,Nprojs);
    imageidx=k;
    im=instack.getImage(k);
    im=double(im);
    
    n=size(im,1);
    if n~=size(im,2)
        error('Images must be square');
    end
    
    [voltage,DefocusU,DefocusV,DefocusAngle,Cs,pixA,A]=...
        cryo_parse_Relion_CTF_struct(CTFdata.data{k});
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
end
outstack.close;
