% Test center of mass correction for noisy projections.
% The script estimates the center of mass for noisy projections at various
% levels of noise and compares the results to the ground truth. This shows
% the senstivity of the estimatation of the center of mass.
% 
% The center of mass is estimated from the noisy projections with and
% without radial masking of the projections.
%
% The projections are centered (not shifts).
% The script generates a table of the estimation error vs SNR. 
%
% Yoel Shkolsniky, November 2014.

clear;
initstate;

K=100;     % Number of projections.
n=65;     % Size of each projection is nxn.
SNRs=[1000 1 1/2 1/4 1/8 1/16 1/32]; % SNRs to test

% Check if using the following masking in mask_projections
mask_radius=floor(n*0.45);
mask_risetime=floor(n*0.1);
center=(n+1)/2;
mask = fuzzymask([n n],2,mask_radius,mask_risetime,[center center]);
% norm(imag(cfft2(mask))) should be tiny

cmref=zeros(K,numel(SNRs),2); % Reference center of mass of each projections.
cmnoisy=zeros(K,numel(SNRs),2); % Center of mass estimated from noisy projections.
cmmasked=zeros(K,numel(SNRs),2); % Center of mass estimated from noisy masked projections.
cmfilt=zeros(K,numel(SNRs),2); % Center of mass estimated from filtered masked projections.

DETAILED=0; % Print estimation results for each image.

if DETAILED
% fprintf('SNR \t\t k \t Ref \t\t Noisy \t\t Noisy err \t masked \t masked err\n');
    fprintf('SNR \t\t k \t noisy \t\t masked  \t filt \n');
end
        
for j=1:numel(SNRs)
    snr=SNRs(j);
    [projs,noisy_projs,~,~]=cryo_gen_projections(n,K,snr,3.1);  % Generate clean projections.
        
    for k=1:K
        masked_proj=noisy_projs(:,:,k).*mask;
        
        clf;
        subplot(1,4,1); imagesc(projs(:,:,k)); colormap(gray); axis image; axis off;
        subplot(1,4,2); imagesc(noisy_projs(:,:,k)); colormap(gray); axis image; axis off;
        subplot(1,4,3); imagesc(masked_proj); colormap(gray); axis image; axis off;
        subplot(1,4,4); imagesc(GaussFilt2(masked_proj,0.3)); colormap(gray); axis image; axis off;
        
        cmref(k,j,:)=CenterOfMass(projs(:,:,k));
        cmnoisy(k,j,:)=CenterOfMass(noisy_projs(:,:,k));        
        cmmasked(k,j,:)=CenterOfMass(masked_proj);
        cmfilt(k,j,:)=CenterOfMass(GaussFilt2(masked_proj,0.3));
        
        if DETAILED
%         fprintf('%6.4e \t %d \t [%+5.2f %+5.2f] \t [%+5.2f %+5.2f] \t %5.3f \t \t[%+5.2f %+5.2f] \t %5.3f\n',...
%             snr,k,cmref(k,j,1),cmref(k,j,2),...
%             cmnoisy(k,j,1),cmnoisy(k,j,2),norm(squeeze(cmref(k,j,:)-cmnoisy(k,j,:))),...
%             cmmasked(k,j,1),cmmasked(k,j,2),norm(squeeze(cmref(k,j,:)-cmmasked(k,j,:))));

            fprintf('%6.4e \t %d \t %6.3f \t %6.3f \t %6.3f\n',...
                snr,k,norm(squeeze(cmref(k,j,:)-cmnoisy(k,j,:))),...
                norm(squeeze(cmref(k,j,:)-cmmasked(k,j,:))),...
                norm(squeeze(cmref(k,j,:)-cmfilt(k,j,:))));
        end
    end
    
end

% Print statistics
%fprintf('SNR \t \t mean noisy \t std noisy \t mean mask \t std mask \t mean filt \t std filt\n');
fprintf('SNR \t \t mean noisy \t mean mask \t mean filt \n');

figure(1); clf;
figure(2); clf;

for j=1:numel(SNRs)
    errnoisy=zeros(K,1);
    errmask=zeros(K,1);
    errfilt=zeros(K,1);
    
    for k=1:K
        errnoisy(k)=norm(squeeze(cmref(k,j,:)-cmnoisy(k,j,:)));
        errmask(k)=norm(squeeze(cmref(k,j,:)-cmmasked(k,j,:)));
        errfilt(k)=norm(squeeze(cmref(k,j,:)-cmfilt(k,j,:)));
    end
    
%     fprintf('%6.4e \t %6.3f \t %6.3f \t %6.3f \t %6.3f \t %6.3f \t %6.3f\n',...
%         SNRs(j),mean(errnoisy),std(errnoisy),mean(errmask),std(errmask),...
%         mean(errfilt),std(errfilt));        
    fprintf('%6.4e \t %6.3f \t %6.3f \t %6.3f\n',...
        SNRs(j),median(errnoisy),median(errmask),median(errfilt));        
    
    figure(1);
    subplot(3,3,j);
    scatter(squeeze(cmnoisy(:,j,1)),squeeze(cmnoisy(:,j,2)));
    axis equal;

    figure(2);
    subplot(3,3,j);
    scatter(squeeze(cmmasked(:,j,1)),squeeze(cmmasked(:,j,2)));
    axis equal;

end
