% test_cryo_align_densities_4
%
% Align two clean density maps where one is rotate/translated relative to
% the other, in the presence of noise. Compare the results to the old code.
%
% Yoel Shkolnisky, July 2016.

initstate;

%% Generate two density maps.
clear;
load cleanrib
vol=real(volref);
vol=GaussFilt(vol,0.8);

[R,~,~]=svd(rand(3)); % Generate random rotation
volRotated=fastrotate3d(vol,R); % Rotate the reference volume by the random rotation
volRotated=reshift_vol(volRotated,(rand(3,1)-1/2)*8); % Shift +/- 4 pixels.

vol=vol/norm(vol(:));
volRotated=volRotated/norm(volRotated(:));

h=figure;
h_old=figure;

SNRs=[100,1,1/2,1/4,1/8,1/16 1/32 1/64 2/128].';
corrsNoisy=zeros(numel(SNRs),1);
corrsClean=zeros(numel(SNRs),1);
corrsNoisy_old=zeros(numel(SNRs),1);
corrsClean_old=zeros(numel(SNRs),1);
resA=zeros(numel(SNRs),1);
resA_old=zeros(numel(SNRs),1);

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
    %% Align using new code
    fprintf('New code SNR=%7.4f:\n',snr);
    verbose=1;

    [Rest,estdx,volAlignedNoisy]=cryo_align_densities(volNoisy,volRotatedNoisy,0,verbose,R,0,200);

    % Correlation of noisy volumes after alignment
    c1=corr(vol(:),volAlignedNoisy(:));
    fprintf('\t Correlation between original and aligned noisy volumes %7.4f\n',c1);
    corrsNoisy(snridx)=c1;
    
    % Correlation of clean volumes according to estimated alignment
    % parameters
    
    volAligned=fastrotate3d(volRotated,Rest.');
    volAligned=reshift_vol(volAligned,estdx);
    
    c2=corr(vol(:),volAligned(:));
    fprintf('\t Correlation between original and aligned clean volumes %7.4f\n',c2);
    corrsClean(snridx)=c2;
    
    % The original and aligned volume should have high agreement according to
    % the FSC. FSC may degrade due to sub-pixel misalignment errors.
    figure(h);
    subplot(3,numel(SNRs),snridx+numel(SNRs));
    plotFSC(vol,volAlignedNoisy,0.5,1,h);
    subplot(3,numel(SNRs),snridx+2*numel(SNRs));
    plotFSC(vol,volAligned,0.5,1,h);
    
    fsc1=FSCorr(vol,volAlignedNoisy);
    resA(snridx)=fscres(fsc1,0.5);
    fprintf('\t Resolution %7.4f\n',resA(snridx));

    %% Align using old code
    fprintf('Old code SNR=%7.4f:\n',snr);
    verbose=1;

    [Rest,estdx,volAlignedNoisy]=cryo_align_densities_old(volNoisy,volRotatedNoisy,0,verbose,R);

    % Correlation of noisy volumes after alignment
    c1=corr(vol(:),volAlignedNoisy(:));
    fprintf('\t Correlation between original and aligned noisy volumes %7.4f\n',c1);
    corrsNoisy_old(snridx)=c1;
    
    % Correlation of clean volumes according to estimated alignment
    % parameters
    
    volAligned=fastrotate3d(volRotated,Rest.');
    volAligned=reshift_vol(volAligned,estdx);
    
    c2=corr(vol(:),volAligned(:));
    fprintf('\t Correlation between original and aligned clean volumes %7.4f\n',c2);
    corrsClean_old(snridx)=c2;
    
    % The original and aligned volume should have high agreement according to
    % the FSC. FSC may degrade due to sub-pixel misalignment errors.
    figure(h_old);
    subplot(3,numel(SNRs),snridx+numel(SNRs));
    plotFSC(vol,volAlignedNoisy,0.5,1,h);
    subplot(3,numel(SNRs),snridx+2*numel(SNRs));
    plotFSC(vol,volAligned,0.5,1,h);
    
    fsc2=FSCorr(vol,volAlignedNoisy);
    resA_old(snridx)=fscres(fsc2,0.5);
    fprintf('\t Resolution %7.4f\n',resA_old(snridx));
    
    
end
fprintf('1/SNR  \t corrNoisy  \t corrNoisy_old  \t corrClean  \t corrClean_old   resA  \t \t resAold\n');
for k=1:numel(SNRs)
    fprintf('%5.0f \t %7.4f  \t %7.4f  \t \t %7.4f   \t %7.4f   \t %7.4f   \t %7.4f\n',...
    1./SNRs(k),corrsNoisy(k), corrsNoisy_old(k), corrsClean(k), corrsClean_old(k), resA(k), resA_old(k));
end