% test_cryo_align_densities_3
%
% Align two clean density maps where one is rotate/translated relative to
% the other, in the presence of noise.
%
% Yoel Shkolnisky, January 2015.
% Revised: July 2016.

initstate;

%% Generate two density maps.

load cleanrib
vol=real(volref);
vol=GaussFilt(vol,0.8);

[R,~,~]=svd(rand(3)); % Generate random rotation
volRotated=fastrotate3d(vol,R); % Rotate the reference volume by the random rotation
volRotated=reshift_vol(volRotated,(rand(3,1)-1/2)*8); % Shift +/- 4 pixels.

vol=vol/norm(vol(:));
volRotated=volRotated/norm(volRotated(:));

h=figure;

SNRs=[100,1,1/2,1/4,1/8,1/16 1/32 1/64 1/128 1/256].';
corrsNoisy=zeros(numel(SNRs),1);
corrsClean=zeros(numel(SNRs),1);

for snridx=1:numel(SNRs)
    snr=SNRs(snridx);
    sigma=sqrt(1/snr);
    n1=sigma*randn(size(vol))./sqrt(numel(vol(:)));
    n2=sigma*randn(size(volRotated))./sqrt(numel(volRotated(:)));
    volNoisy=vol+n1;
    volRotatedNoisy=volRotated+n2;
    
        
    %% Visulaize noisy level
    figure(h);
    subplot(3,numel(SNRs),snridx);
    imagesc(volNoisy(:,:,floor((size(volNoisy,3)+1)/2))); 
    colormap(gray);
    axis image;
    axis off;
    %% Align
    verbose=1;

    [Rest,estdx,volAlignedNoisy]=cryo_align_densities(volNoisy,volRotatedNoisy,0,verbose,R,0,50);

    % Correlation of noisy volumes after alignment
    c1=corr(vol(:),volAlignedNoisy(:));
    fprintf('Correlation between original and aligned noisy volumes %7.4f\n',c1);
    corrsNoisy(snridx)=c1;
    
    % Correlation of clean volumes according to estimated alignment
    % parameters
    
    volAligned=fastrotate3d(volRotated,Rest.');
    volAligned=reshift_vol(volAligned,estdx);
    
    c2=corr(vol(:),volAligned(:));
    fprintf('Correlation between original and aligned clean volumes %7.4f\n',c2);
    corrsClean(snridx)=c2;
    
    % The original and aligned volume should have high agreement according to
    % the FSC. FSC may degrade due to sub-pixel misalignment errors.
    figure(h);
    subplot(3,numel(SNRs),snridx+numel(SNRs));
    plotFSC(vol,volAlignedNoisy,0.5,1,h);
    subplot(3,numel(SNRs),snridx+2*numel(SNRs));
    plotFSC(vol,volAligned,0.5,1,h);
    
end
fprintf('     1/SNR  corrNoisy  corrClean\n');
disp([1./SNRs corrsNoisy corrsClean]);