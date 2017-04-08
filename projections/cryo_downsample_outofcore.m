function cryo_downsample_outofcore(instackname,outstackname,szout,verbose)
% CRYO_DOWNSAMPLE_OUTOFCORE      Downsample projections
%
% cryo_downsample_outofcore(instackname,outstackname,szout)
%   Downsample projections in the MRC file named instackname, and write the
%   downsampled images to the MRC file outstackname. 
%   The images are downsampled to size szout, which is either a scalar or a
%   vector of the dimension of the output images.
%   The function does not load the entire stack into memory.
%
% Example:
%   cryo_downsample_outofcore('instack.mrc','outstack.mrc',[65 65]);
%
% Yoel Shkolnisky, May 2016.

if nargin<4
    verbose=1;
end

instack=imagestackReader(instackname,100);
Nprojs=instack.dim(3);
outstack=imagestackWriter(outstackname,Nprojs,1,100);

if verbose
    printProgressBarHeader;
end
for k=1:Nprojs
    if verbose
        progressTicFor(k,Nprojs);
    end
    im=instack.getImage(k);    
    downsampled_im=cryo_downsample(im,szout,0);  % We are downsampling a 
        % single image so set the stack flag to 0.
    downsampled_im=single(downsampled_im);
    outstack.append(downsampled_im);
end
outstack.close;