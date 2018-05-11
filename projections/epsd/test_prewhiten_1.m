
N=129; % Each noise image is of size NxN.
Ks=1000; % Number of noise images to generate.

initstate;

% Generate a stack of noise images
[noise,Sfft,Sref,T1]=noise_exp2d(N,max(Ks),1);

max_d=floor(N/3);
[P2,~,~,~]=cryo_epsdS(noise(:,:,1:Ks),1:N^2,max_d,1);

% Display power spectrum
figure;
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
title(sprintf('Central profile of power spectra (K=%d)',Ks));


% Prewhiten
prewhitened=cryo_prewhiten(noise, P2);
[prewhitenedP2,~,~,~]=cryo_epsdS(prewhitened,1:N^2,max_d,1);

subplot(2,2,2)
prewhitenedP1d=prewhitenedP2(N,:);    % Central 1D profile of PSD estimated using cryo_epsdS.
plot(prewhitenedP1d,'x-','Color','b');
ylim([0, 1.1]);
title('Prewhitened spectrum');

figunits=get(gcf,'units');
figsize=get(gcf,'outerposition');
set(gcf,'units','normalized','outerposition',[0 0 1 1])

figname=sprintf('prewhiten_%d_%d',N,Ks);
hgsave(figname);

set(gcf,'units',figunits);
set(gcf,'outerposition',figsize);
