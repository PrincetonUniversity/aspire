function [vol,projs]=cryo_orient_projections_auxpreprocess(vol,projs)
%
% Preprocess the volume and projections for the various
% cryo_orient_projections_* functions.
%
% The preprocessing consists of:
% 1. Downsample vol to the dimensions of projs.
% 2. Estimate the SNR of the projections. If the SNR is low (below 100),
%    estimate the power spectrum of the noise in the projections,prewhiten
%    them, and normalize them to background mean 0 and variance 1.
%
% Yoel Shkolnisky, August 2015.

szvol=size(vol);

% Resample reference volume to the size of the projections
szprojs=size(projs); % Check that projections are square.
if szprojs(1)~=szprojs(2)
    error('Projections to orient must be square.');
end

if any(szvol-szvol(1)) % Check that the volume has all dimensions equal.
    error('Volume must have all dimensions equal');
end

% Mask and filter the volume
n=size(vol,1);
vol=vol.*fuzzymask(n,3,floor(0.45*n),floor(0.05*n));
%vol=GaussFilt(vol,0.3);

% Mask and filter the projections
for k=1:size(projs,3)
    p=projs(:,:,k);
    p=p.*fuzzymask(n,2,floor(0.45*n),floor(0.05*n));
    %p=GaussFilt(p,0.3);
    projs(:,:,k)=p;
end

projs=cryo_globalphaseflip(projs);

% n=szprojs(1);
% log_message('Resampling refence volume from %d to %d',szvol(1),n);
% vol=cryo_downsample(vol,[n,n,n],0);
% 
% snr=cryo_estimate_snr(projs);
% log_message('Estimated SNR of projections = %5.2f',snr);
% if snr<100 % Don't bother when the SNR is high
%     % Prewhiten projections
%     log_message('Estimating noise power spectrum');
%     n=size(projs,1);
%     log_message('Each projection of size %dx%d',n,n);
%     psd = cryo_noise_estimation(projs);
%     log_message('Finished noise power spectrum estimation');
% 
%     log_message('Prewhitening images');
%     projs = cryo_prewhiten(projs, psd);
%     log_message('Finished prewhitening images');
% 
%     % Normalize projections background
%     projs=cryo_normalize_background(projs,floor(szvol(1)/2),0);
% end

% % Compute frequnecy support of the volume.
% ravg3d=cryo_radial_powerspect_3d(vol);
% cs=cumsum(ravg3d);
% th=1-10^(-5);
% idx=find(cs>th,1,'first'); % The cutoff frequnecy is set to be the index 
%                              % where the cumulative enery is above th.
% cutoff3d=idx/numel(ravg3d) ; % cutoff frequnecy in normalized units.
% log_message('Cutoff frequnecy of the volume %5.2f (cumulative energy above %5.2f)',cutoff3d,th);
% 
% % Compute frequnecy support of the projections.
% ravg2d=cryo_radial_powerspect_2d(projs);
% cs=cumsum(ravg2d);
% th=1-10^(-5);
% idx=find(cs>th,1,'first'); % The cutoff frequnecy is set to be the index 
%                              % where the cumulative enery is above th.
% cutoff2d=idx/numel(ravg2d) ; % cutoff frequnecy in normalized units.
% log_message('Cutoff frequnecy of the projections %5.2f (cumulative energy above %5.2f)',cutoff2d,th);
% 
% % Filter both projections and volume
% cutoff=min(cutoff2d,cutoff3d);
% log_message('Filtering volume and projections to cutoff %5.2f',cutoff);
% vol=GaussFilt(vol,cutoff);
% projs=GaussFilt(projs,cutoff,1); % 1 means we have a stack of images.