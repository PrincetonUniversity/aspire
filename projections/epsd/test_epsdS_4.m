% Test the function cryo_epsdS.
% 
% The script runs the following steps:
% 1. Generate colored-noise images with known power spectrum.
% 2. Estimate the isotropic power spectrum of the noise using cryo_epsdS.
% 3. Compare the known and estimated spectra.
%
% The code works for both even and odd images sizes N.
%
% Yoe Shklnisky, October 2014.

N=129; % Each noise image is of size NxN.
Ks=[100 1000 10000]; % Number of noise images to generate.

initstate;

% Generate a stack of noise images
[noise,Sfft,Sref]=noise_fakectf2(N,max(Ks));

for j=1:numel(Ks)
    fprintf('Processing K=%d\n',Ks(j));
    
    % Estimate power spectrum of noise images    
    max_d=floor(N/2);
    [P2,R,~,x]=cryo_epsdS(noise(:,:,1:Ks(j)),1:N^2,max_d,1);
    
    if j==1 % Allocate memory for all power spectra of all tests.
        xs=zeros(numel(R),numel(Ks));
        Rs=zeros(numel(R),numel(Ks));
        P2s=zeros(2*N-1,2*N-1,numel(Ks));
    end

    P2 = P2/norm(P2(:));
    
    Rs(:,j)=R;
    P2s(:,:,j)=P2;
    
end

% Display power spectrum
for j=1:numel(Ks)
    figure;
    P2=P2s(:,:,j);      % 2D power spectrum estimated using cryo_epsdS.
    P1dSref=Sref(N,:);  % Central 1D profile of the analytic power spectrum. 
    P1dSfft=Sfft(N,:);  % Central 1D profile estimated using FFT on noise-only images.
    P1dSest=P2(N,:);    % Central 1D profile of PSD estimated using cryo_epsdS.

    subplot(1,2,1); % Plot the profiles of the 3 PSDs.
    hold on;
    plot(P1dSref,'.-','Color','g');
    plot(P1dSfft,'-','Color','m','LineWidth',2);
    plot(P1dSest,'x-','Color','b');
    legend('Sref',...
        sprintf('Sfft (%5.3f)',corr(Sfft(:),Sref(:))),...
        sprintf('cryo_epsdS  (%5.3f)',corr(P2(:),Sref(:))));
    xlabel('frequency','FontSize',12);
    ylabel('power spectrum','Rotation',90,'FontSize',12);
    set(gca,'FontSize',12);
    hold off;
    title(sprintf('Central profile of power spectra (K=%d)',Ks(j)));
    
    subplot(1,2,2) % Relative error between central 1D profile analytic and 
                   % estimated (using cryo_epsdS) power spectram
    plot(P1dSref-P1dSest);
    title(sprintf('Relative L2 err=%6.4f',norm(P2(:)-Sref(:))/norm(Sref(:))));
    
    figunits=get(gcf,'units');
    figsize=get(gcf,'outerposition');
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    
    figname=sprintf('fakectf2_%d_%d',N,Ks(j));
    hgsave(figname);
    
    set(gcf,'units',figunits);
    set(gcf,'outerposition',figsize);
end
