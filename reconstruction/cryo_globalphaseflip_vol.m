function [vol,doflip,signalmean,noisemean]=cryo_globalphaseflip_vol(vol)
%
% CRYO_GLOBALPHASEFLIP_VOL  Apply global phaseflip to a volume
%
% [vol,doflip,signalmean,noisemean]=cryo_globalphaseflip(vol)
%   Check if a volume should be phase flipped so that the volume
%   corresponds brigther pixels and the background corresponds 
%   to darker pixels. This is done by comapring the mean in a small circle
%   around the origin (supposed to corrspond to the volume) with the mean
%   of the noise, and making sure that the mean of the volume is larger.
%
%   Input parameters:
%       vol  3D volume. 
%
%   Output parameters:
%       vol         (flipped) volume.
%       doflip      1 if phase was flipped and 0 otherwise.
%       signalmean  Mean of the signal.
%       noisemean   Mean of the noise.
%
%   Examples:
%       stack=ReadMrc('vol.mrc');
%       flippedvol=cryo_globalphaseflip_vol(vol);
%
% Yoel Shkolnisky, April 2017

sz=size(vol);

if numel(sz)~=3
    error('First input must be a volme (3D array)');
end

if sz(1)~=sz(2) || sz(1)~=sz(3)
    error('Volume must have all dimensions equal');
end

n=sz(1);
center=(n+1)/2;
[I,J,K]=meshgrid(1:n,1:n,1:n);
r=sqrt((I-center).^2+(J-center).^2+(K-center).^2);
sigind=r<round(n/4); % Indices of signal samples
noiseind=r>round(n/2*0.8);

signalmean=mean(vol(sigind));
noisemean=mean(vol(noiseind));        

doflip=0;
if signalmean<noisemean
        doflip=1;
        vol=-vol;
end
