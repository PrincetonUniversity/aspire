
function [res,fc] = fsc( v1, v2, do_mask, pixelsize )
%FSC Mask density maps and compute the Fourier Shell Correlation between them
%
% This is a wrapper of FSCorr, written by Ido Greenberg, 2016.

if ~exist('do_mask','var'); do_mask=true; end

n = size(v1, 1);
if ~exist('pixelsize','var')
    % Default pixel size correponds to current 80s data
    %pixelsize = 1.5*247/n;
    pixelsize = 1.34*359/n;
end

% masking
if do_mask
    rmask1 = cryo_masking_radius_3d(v1,0.99,0);
    rmask2 = cryo_masking_radius_3d(v2,0.99,0);
    rmask = min( max(rmask1,rmask2) , (n+1)/2 );
    log_message('Masking radius [pixels]: %d/%d', rmask, (n+1)/2);
    v1 = cryo_mask_volume(v1,rmask);
    v2 = cryo_mask_volume(v2,rmask);
end

% FSC
fc = FSCorr(v1,v2);

res_f = fscres(fc,0.143);
res = 2*pixelsize*numel(fc)/res_f;

end
