function WriteJpeg(m,filename,clipThreshold)
% function WriteJpeg(m,filename,clipThreshold)
% Autoscale the image m and write it out as a jpeg image.  If filename has
% no extension add '.jpg' to it.  clipThreshold is the fraction of gray
% values that are clipped at black and white ends of the range; default is
% 1e-3.
if nargin<3
    clipThreshold=1e-3;
end;
if numel(regexp(filename,'.+\.jpg'))+numel(regexp(filename,'.+\.jpeg'))==0
    filename=[filename '.jpg'];
end;
imwrite(uint8(imscale(rot90(m),256,clipThreshold)),filename);
