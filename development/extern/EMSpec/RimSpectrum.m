function sp=RimSpectrum(m0,msk,displayOn)
% function sp=RimSpectrum(m0,msk)
% Compute the power spectrum from a stack of images m0. The mask msk is
% used to mask out the central part of each image, where the particle is.
% For example, to plot the spectrum of a stack of particle images,
%   n=size(m0,1);
%   msk=1-fuzzymask(n,2,n/2-4,2);  % rim thickness is 4 pixels
%   sp=RimSpectrum(m0,msk);
%    % The returned power spectrum has n/2 points.
%    % pixA is the angstroms per pixel
%   semilogy((0:n/2-1)/(n*pixA),sp);
% Scaling is such that an image with unity variance white pixel noise gives
% a flat spectrum (except for the first, zero-frequency point) of value 1.
%
% Algorithm: from masked images, we compute the power spectrum, and do ifft
% and radial averaging to get an autocorrelation function.  We also compute
% the acf of the mask.  Then we divide these two acfs and do a Hankel
% transform to get back to a radially-averaged power spectrum.

testMode=0;
    sr=2; sc=3;  % subplot size

if testMode
    %% make fake data
    n=128;  % Basic image size
    n2=2*n; % Padded image size
    ni=2048; % Number of images
    
    mw=2;  % mask transition width
    msk=8;  % mask thickness
%     nargin=3;
    displayOn=1;
    
    % % base=fftshift(fuzzymask(n2,2,0.45*n2,.05*n2) .*30./((Radius(n2))+10));
    base=fftshift(fuzzymask(n2,2,0.45*n2,.05*n2) .*10./(sqrt(Radius(n2))+3));
    % % H=fftshift(sqrt(Radius(n2))/n.*CTF(n2,2,.025,.5,2,50,.2));
    H=fftshift(CTF(n2,2,.025,.5,2,50,.2));
    figure(1); SetGrayscale; subplot(sr,sc,1);
    imacs(fftshift(H));
    subplot(sr,sc,2);
    imacs(fftshift(sqrt(base)));
    %
    mt=msk;
    msk=Crop(1-fuzzymask(n,2,n/2-mt,mw),n2);
    %  msk=Crop(1-fuzzymask(n,2,0,4),n2);
    subplot(sr,sc,3); imacs(msk);
    drawnow;
    m0=zeros(n,n,ni);
    
    disp('Making images');
    if testMode>1
        for i=1:ni
            %     m1=Crop(randn(n/2,n/2),n);  % padded random image
            mf=fftn(randn(n2,n2));  % We'll mask out the center of a double-sized image
            mf2=fftn(randn(n2,n2));
            %     m(:,:,i)=real(ifftn(mf.*H+mf2.*base)).*msk;
            m0(:,:,i)=Crop(real(ifftn(mf.*H+mf2.*base)),n);
        end;
    else
        for i=1:ni
            %       m0(:,:,i)=GaussFilt(randn(n,n),.5);
            m0(:,:,i)=randn(n,n);
        end;
        
    end;
end;
%% Actual code starts here
if nargin<3
    displayOn=0;
end;
n=size(m0,1);
ni=size(m0,3);
n2=2*n;

if numel(msk)==1  % it's a scalar, sets the rim width
    mt=max(msk/3,1);  % rim risetime
    msk=1-fuzzymask(n,2,n/2-msk,mt);  % create a mask
end;
msk=Crop(msk,n2);  % Pad the mask to double size

% % Computing full power spectrum
% rsp=zeros(n/2,1);
% for i=1:ni
%     rsp=rsp+RadialPowerSpectrum(m0(:,:,i));
% end;
% rsp=rsp/ni;  % should be unity for white noise
%

% Make the masked images
m=zeros(n2,n2,ni);
smsk=sum(msk(:));
for i=1:ni
    m1=Crop(m0(:,:,i),n2).*msk;  % pad the entire stack
    m(:,:,i)=m1-msk.*(sum(m1(:))/smsk); % Force zero mean
end;

if displayOn
    figure(1);
    SetGrayscale;
    NumberOfImages=ni
    subplot(sr,sc,2); imacs(m(:,:,1));
    drawnow;
    disp('Computing spectra');
end;
%%
% accumulate the spectra
% sp=zeros(n2,n2);
sps=zeros(n2,n2,ni);
for i=1:ni
    %     sp=sp+abs(fftn(m(:,:,i))).^2;
    sps(:,:,i)=abs(fftn(m(:,:,i))).^2;
end;
% sp=sp/(ni*n^2);
sp=median(sps,3)/n^2;  % less sensitive to outliers.
% sp=mean(sps,3)/(n^2);

spm=abs(fftn(msk)).^2/n^2;  % Spectrum of mask function
if displayOn
    subplot(sr,sc,4);
    imacs(fftshift(sp).^.2);
    %     subplot(sr,sc,6);
    %     plot([Radial(fftshift(sp)) 1e7*Radial(fftshift(H.^2))]);
    subplot(sr,sc,5);
    imacs(fftshift(spm).^.2);
    drawnow;
    %%
    disp('Computing autocorrelations');
end;
% autocorrelations
% padded spectrum by os
os=2;  % os=4 is nice, os=2 gives essentially the same results.
spx=fftshift(Crop(fftshift(sp),n2*os));  % pad the spectrum
c=real(ifftn(spx));                      % get the raw acf
spmx=fftshift(Crop(fftshift(spm),n2*os));  % pad the mask spectrum
cm=real(ifftn(spmx));                      % acf of mask.

if displayOn
    subplot(sr,sc,4);
    imacs(fftshift(abs(c)).^.2);
    subplot(sr,sc,5);
    imacs(fftshift((cm)));
    drawnow;
end;
% cr=sectr(fftshift(c));

% cr=Radial(fftshift(c.*cm));
% cmr=Radial(fftshift(cm.^2));
% The gridding-based polar conversion works much better than the
% nearest-neighbor code in Radial().
cr=mean(gridToPolar(fftshift(c)),2);
cmr=mean(gridToPolar(fftshift(cm)),2);
crr=cr./cmr;  % Corrected autocorrelation
if displayOn
    subplot(sr,sc,1);
    disfac=cmr(1)/cr(1);
    q=[cr cmr/disfac cmr(1)*crr];
    % plot(q(1:n*os/4,:));
    plot(crr);
    disp('Hankel transform');
end;

% Hankel
na=numel(crr);
accum=zeros(na,1);
rp=(0:na-1)'*pi/na;  % Radial position values, times frequency scaler
for i=1:na
    s=(i-1);  % radial frequency value
    accum=accum+besselj(0,rp*s)*s*crr(i);
end;
s0=accum*8/(os^2*1.08)*sqrt(2);  % sqrt(2) for median spectrum calculation.
    % I don't know where the 1.05 arises.


% %
% % % Brute-force code
% %
% % %
% % csym=ToRect(crr);
% csym=gridToRect(repmat(crr,1,numel(crr)*6));  % symmetrized acf
% subplot(sr,sc,5);
% imacs(Crop(csym,n/4*os));
% % %
% % % h0=zeros(n*os,1);
% % h0=sectr(fftshift(H.^2));
% ssym=real(fftn(fftshift(csym)));
% sb0=Radial(fftshift(ssym));
% subplot(sr,sc,6);
% plot(sb0(1:round(n*.9)));
%

%  h0=Radial(fftshift(H.^2+base.^2));

% s0=GaussFilt(s0,.1)/n/4;
ss0=max(s0,0);
ss0=.5*(ss0(2:2:n-1)+ss0(3:2:n));
ss0(n/2)=ss0(n/2-1);
ss0=max(ss0,0);
% q=[ss0(1:n) n*os*sqrt(os/2)*h0];
% q=[ss0(1:n)];
sp=ss0;

if displayOn
    subplot(sr,sc,3);
    np=round(.9*n/2);
    np=n/2;
    plot((0:np-1)/n,sp(1:np,:),'k.-');
    subplot(sr,sc,6);
    np=floor(.45*n);
    semilogy((0:np-1)/n,sp(1:np,:));
    s0(1:8)
end;
% %%
%
%
%
%
% subplot(sr,sc,4);
% plot(accum);
% subplot(sr,sc,3);
% plot([2e4*H(1:128,n+1).^2 1e0*accum]);
