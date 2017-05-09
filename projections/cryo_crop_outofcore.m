function cryo_crop_outofcore(instackname,outstackname,croppeddim)
% CRYO_CROP_OUTOFCORE      Crop projections
%
% cryo_crop_outofcore(instackname,outstackname,croppeddim)
%   Crop to the projections in the MRC file instackname to size
%   croppeddim. The cropped images are written to the MRC file
%   outstackname.
%   croppeddim is either a scalar or a vector of the dimension of the
%   output images.
%
% Example:
%   cryo_phaseflip_outofcore('instack.mrc','outstack.mrc',65); 
%   cryo_phaseflip_outofcore('instack.mrc','outstack.mrc',[65 65]); 
%
% Yoel Shkolnisky, May 2016.

instack=imagestackReader(instackname,100);
Nprojs=instack.dim(3);
outstack=imagestackWriter(outstackname,Nprojs,1,100);

%printProgressBarHeader;
for k=1:Nprojs
%    progressTic(k,Nprojs);
    im=instack.getImage(k);
    cropped_im=cryo_crop(im,croppeddim,1); 
    cropped_im=single(cropped_im);
    outstack.append(cropped_im);
end
outstack.close;
