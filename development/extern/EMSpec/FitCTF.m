function [P c spect]=FitCTF(mc,Pa,pixA,maxres,minres,test,options)
% function [P c spect]=FitCTF(mc,Pa,pixA,maxres,minres,test)
%  Derived from ppcfit2, uses the algorithm of N. Grigorieff's CTFFIND3.
% Determine the CTF parameters P from the image mc.  The cross-correlation
% coefficient c is returned too, as is the averaged spectrum spect (zero
% frequency in the center).
% Positive defocus values correspond to *underfocus*.
% The image is assumed to be pre-whitened.
% Parameter ranges are defined by the structure Pa.  If vectors of values
% are given, all combinations of these values are tested in the brute-force
% search.  Afterwards, parameters are optimized using Simplex.
% The arguments pixA, maxres, minres are all in angstroms.
% The test argument controls display of results.  test=1 means graphical
% output (default).  test=2 gives text display of execution times.
% options is a struct of parameters for special options with fields
% spFromWholeMicrograph: use entire micrograph for spectrum.
%
% Here is an example:
%   Pa.lambda=EWavelength(200); % lambda in angstroms
%   Pa.defocus=0:0.5:20;  % in microns
%   Pa.deltadef=-2:.4:2;  % astigmatism
%   Pa.theta=0;
%   Pa.alpha=.07;  % This can also be fitted, for use with phase plate.
%       It is actually an angle in radians; to fit it, give a vector of
%       values e.g. 0.5:0.5:2.
%    % The following parameters are not fitted
%   Pa.Cs=2;
%   Pa.B=200;  % This is not fitted
%    % Actual function call
%   P=ctfit2(m,Pa,1.4,8,50);
%   P      % print out the values
% Tweaked to work better with Liguo's data -fs 20 Apr 12
%  - assumes flat HF spectrum (CCD prewhitening)
%  - no hp filtering (drift removal) of spectrum
%  - multiply spectrum by f

if nargin<7
    options=struct;
end;
if nargin<6
    test=1;  % by default, show graphics.
end;

% Set up default options
if ~isfield(options,'spFromWholeMicrograph')
    options.spFromWholeMicrograph=0;
end;

NumRefs=2;
persistent R;  % We cache the references to allow faster execution of repeats.

mc=single(mc);
mc=mc-mean(mc(:));
[nx ny]=size(mc);
%
% maxres=6;
% minres=40;

Pa.pixA=pixA;  % pixel size in A.

if nx>5000
    nu=256;  % FFT block size, must be a multiple of 4.
else
    nu=256;
end;

fExponent=1;

% code for hp filtering of spectrum, presently unused.
defocusMin=min(Pa.defocus);
firstZero=sqrt(defocusMin*Pa.lambda*1e4); % first zero, in Å
w0=firstZero/2;  % High-pass filter for spectrum, in A.  Should be smaller than
w0=100;
w=w0/pixA;     % width of real-space kernel.
% w=0;
disexp=0.2; % display exponent
% pwExp1=4;  % pre-whitening positive exponent
% pwExp4=2; % pre-whitening negative exponent (r^4)


df=1/(pixA*nu);  % spatial frequency unit
freqs=(0:df:df*(nu/2-1))';

% Compute local power spectra
[nx ny]=size(mc);
nv=round(nu/4);   % overlap
tx=TileCoords(nx,nu,nv);  % no. of tiles in X direction
ty=TileCoords(ny,nu,nv);  % no. of tiles in Y direction

if test>0  % Show the micrograph
    figure(1);
    SetGrayscale;
    subplot(2,3,1);
    imacs((1:nx)*pixA/10,(1:ny)*pixA/10,mc);
    xlabel('nm');
    drawnow;
end;

sd=zeros(tx,ty);
% Make a window for spectral estiation
window=fuzzymask(nu,2,0.45*nu,0.1*nu);
winsum=sum(window(:));

sps=zeros(nu,nu,tx,ty);

if test>1
    disp('Computing spectra...');
    tic
end;
for ix=1:tx
    [x0 x1 u1]=TileCoords(nx,nu,nv,ix);
    
    for iy=1:ty
        [y0 y1 v1]=TileCoords(ny,nu,nv,iy);
        % I don't think it's necessary to remove gradients, but Niko does
        % this:
        tm=RemoveGradients(double(mc(x0+1:x0+nu,y0+1:y0+nu)));
        tm=tm.*window;
        tm=tm-sum(tm(:))*window/winsum;
        
        % subplot(2,1,1); imacs(tm); title([ix iy]);
        
        sp2=abs(fftn(tm)).^2;
        % subplot(2,1,2); imacs(fftshift(sqrt(sp2)));
        % drawnow;
        sps(:,:,ix,iy)=sp2;
        sd(ix,iy)=sqrt(sum(sp2(:)))/nu;
    end;
end;
if test>1
    toc
end;



sdmn=min(sd(:));
sdmx=max(sd(:));
if options.spFromWholeMicrograph
    sdmin=sdmn;
    sdmax=sdmx;
else
    % Make a histogram of the s.d. values and find tiles having s.d. near the
    % mode.
    nbins=tx+2;
    hthresh=0.2;
    % hthresh=-1;  % no selection
    
    dbin=(sdmx-sdmn)/(nbins-1);
    bins=sdmn:dbin:sdmx;
    
    [h x]=histc(sd(:),bins);
    
    [mxh j]=max(h);
    % j is the mode of the distribution
    % Find the first bin below hthresh on the high side
    jmax=j-1+find(h(j:nbins)<hthresh*mxh,1,'first');
    if numel(jmax)<1
        jmax=nbins;
    end;
    
    % find the bin below hthresh on the low side
    jmin=find(h(1:j)<hthresh*mxh,1,'last');
    if numel(jmin)<1
        jmin=1;
    end;
    sdmin=bins(jmin);
    sdmax=bins(jmax)+dbin;
end;
% Show histogram
% subplot(2,3,2);
% bar(x,h);
% hold on
% plot(sdmin,0,'w.');
% plot(sdmax,0,'w.');
% hold off

% Show regions that we used.
if test
    subplot(2,3,2);
    imacs((sd<=sdmax).*(sd>=sdmin));
    drawnow;
end;

cumsp=zeros(nu,nu);
count=0;
for ix=1:tx
    for iy=1:ty
        if (sd(ix,iy)<=sdmax) && (sd(ix,iy)>=sdmin)
            cumsp=cumsp+sps(:,:,ix,iy);
            count=count+1;
        end;
    end;
end;
if test>1
    disp([num2str(count) ' tiles of ' num2str(tx*ty) ' used']);
end;
spect=fftshift(single(cumsp/count));
% pre-whitening correction
% r=fftshift(Radius(nu))/nu;
% cumsp0=cumsp; %%%
% cumsp=cumsp.*exp(r.*(pwExp1-pwExp4*r.^3));
%
% remove the crystal spots
% cumsp=fftshift(RemoveSpots(fftshift(cumsp),[3.5 3]));

% Filter the square root of the spectrum.
sqrsp=sqrt(cumsp);
if test
    subplot(2,3,4);
    xs=(-nu/2:nu/2-1)*df;
    imacs(xs,xs,fftshift(sqrsp.^(disexp*2)));
    xlabel('A^{-1}');
%     LabelImage(nu,df,'A^{-1}');
end;

% % We want to smooth the spectrum to about firstZero/6
% %  width will be .16/(firstZero*df) so frequency is 1/width
% wl=df*(firstZero/16).^2*80;  % stronger filtering for smaller defocus
%
% H=(1-fuzzymask(nu,2,w,w/5)).*fuzzymask(nu,2,wl,wl/2);  % filter for spectrum
% % convolve with kernel
%
% H=fuzzymask(nu,2,wl*nu,wl*nu/4);
% filtsp=real(ifftn(fftn(sqrsp).*ifftshift(H)));  % gives lots of NaNs.  Why?
%
% diffsp=fftshift(filtsp);

%%%%%%%%%%%%%%%%%%%%
diffsp=fftshift((sqrsp));  %% no filtering at all.
% diffsp=fftshift(cumsp./filtsp-filtsp);
% diffsp=fftshift(cumsp);

% look at the subtracted spectrum

% % subplot(2,3,3);  % Compare the smoothed and original spectra
% % radSpecs=[Radial(fftshift(filtsp))' Radial(fftshift(sqrt(cumsp)))'];
% % semilogy(freqs,radSpecs);
% % xlabel('A^{-1}');
% % drawnow;

radsp=Radial(fftshift(cumsp));


% determine resolution limits

rmin=1/(df*minres);
rmax=1/(df*maxres);
outerlimits=fuzzymask(nu,2,rmax,rmax/10);
limits=outerlimits-fuzzymask(nu,2,rmin,1);
% zero the baseline of diffsp
annulus=outerlimits-fuzzymask(nu,2,rmax*.6,rmax/10);
diffsp=diffsp-(annulus(:)'*diffsp(:))/sum(annulus(:));
diffsp=diffsp.*(Radius(nu)/nu).^fExponent;  % scale up hf part
if test
    subplot(2,3,5);
    xs=(-nu/2:nu/2-1)*df;
    imacs(xs,xs,limits.*diffsp);
    xlabel('A^{-1}');
    
    subplot(2,3,6);
    rdiff=Radial(diffsp.*limits);
    % % plot(freqs,rdiff);
    drawnow;
end;

% Pao=Pa;
% Pao.defocus=-Pa.defocus;

halflimits=limits(nu/2+1:nu,:);
halflimits(1,nu/2:nu)=0;  % zero out the upper meridian

halfsp=diffsp(nu/2+1:nu,:).*halflimits;
halfspv=halfsp(:);
halfspv=halfspv/sqrt(halfspv'*halfspv);

% % Here we test for an existing set of references, and compute new ones
% % only if necessary.
% First, check whether the static variable R exists at all.
b=1;
try
    b=numel(R);
catch
    b=0;
end;

% Now see if the problem has changed from previous calls.
ind=0;
for i=1:b
    if StructsEqual(R(i).Pa0,Pa)
        ind=i;
    end;
end;

if ind==0
    ind=min(b+1,NumRefs);
    if test>1
        disp('Making new CTF references');
    end;
    [R(ind).refs R(ind).refsz]=MakeCTFRefs(nu,pixA,Pa,halflimits);
    %     [R(ind).refso R(ind).refsz]=MakeCTFRefs(nu,res,Pao,halflimits);
    
    R(ind).refs=reshape(R(ind).refs,nu^2/2,prod(R(ind).refsz));
    %     R(ind).refso=reshape(R(ind).refso,nu^2/2,prod(R(ind).refsz));
    R(ind).Pa0=Pa;
    if test>1
        whos R
    end;
end;

if test>1
    disp('Cross-correlation')
end;

cc=halfspv'*R(ind).refs;  % Cross correlation done here!

% pick up the indices of the best match
[mxc mxi]=max(cc);
[m l i j k]=ind2sub(R(ind).refsz,mxi);

if test>1
    Initial_c=mxc
end;

Ps=Pa;
P=Ps;

% Set the starting search values
P.defocus=Ps.defocus(i);
P.deltadef=Ps.deltadef(j);
P.theta=Ps.theta(k);
P.alpha=Ps.alpha(l);
P.B=P.B(m);
% P

if test
    % Make the split-power spectrum display
    subplot(2,3,5);
    dspect=imscale(limits.*diffsp);
    dref=imscale(limits.*CTF(nu,pixA,P).^2);
    dspect(1:nu/2+1,:)=dref(1:nu/2+1,:);
    xs=(-nu/4:nu/4-1)*df;
    imac(xs,xs,Crop(dspect,nu/2));
    xlabel('A^{-1}');
    
    subplot(2,3,4);
    xs=(-nu/2:nu/2-1)*df;
    imac(xs,xs,dspect);
    xlabel('A^{-1}');

    
    drawnow;
end;

ccx=reshape(cc,R(ind).refsz);

% make a contour plot of alpha vs defocus
% figure(3); clf;
% colormap jet
if test
    subplot(2,3,3);
    if numel(Ps.alpha)>1
        contourf(Ps.alpha,Ps.defocus,squeeze(ccx(m,:,:,j,k))',10);
        xlabel('Alpha');
        ylabel('Defocus');
        colorbar;
        title('Underfocus is positive');
    end;
    subplot(2,3,3);
    contourf(Ps.deltadef,Ps.defocus,squeeze(ccx(m,l,:,:,k)),10);
    xlabel('Delta-defocus');
    ylabel('Defocus');
    colorbar;
    drawnow;
end;
% Do the optimization
% disp('Final optimization');
if test>1
    disp('Final optimization.');
end;
P.theta=0;
fitMask=ones(4,1);  % determines which parameters are to be fitted.
if numel(Ps.alpha)<2  % Don't fit alpha
    fitMask(4)=0;
    niters=80;
else
    niters=160;
end;
p=[P.defocus; P.deltadef; P.theta; P.alpha];
alphaStep=max(max(abs(P.alpha))/4,0.03);

for is=1:2  % do extra restarts of the optimization, to get out local optima.
    %     disp(p');
    p=Simplex('init',p,[0.1 0.1 0.1 alphaStep],fitMask);
    for i=1:80
        P.defocus=p(1);
        P.deltadef=p(2);
        P.theta=p(3);
        P.alpha=p(4);
        c=CTFCrossCorr(diffsp,pixA,limits,P,0);
        p=Simplex(-c);
    end;
    % disp([p(1:2)' 180/pi*p(3:4)']);
end;

if test
%     subplot(2,3,4);
%     dspect=imscale(fftshift(sqrsp.^(disexp*2)));
%     dspect=(imscale(ImEqualize(GaussFilt(fftshift(sqrsp),0.2)).^2)*1+0);
%     [dref mr ma]=imscale(limits.*CTF(nu,pixA,P).^2);
%     dref=ma+mr*CTF(nu,pixA,P).^2;
%     % dref=ImEqualize(CTF(nu,res,P).^2);
%     dspect(1:nu/2+1,:)=dref(1:nu/2+1,:);
%     imac(Crop(dspect,nu/2));

%       circularly-averaged spectrum display--not successful.
%     annMask=fuzzymask(nu,2,nu*0.45,nu*.01)-fuzzymask(nu,2,nu*.49,nu*.01);
%     ssqrsp=fftshift(sqrsp);
%     sqrspm=ssqrsp-(annMask(:)'*ssqrsp(:)/sum(annMask(:)));
%     imac(Crop(imscale(abs(RotSymmetrize(sqrspm,4)).^(disexp*2),256,.01),nu/2));
%     LabelImage(nu/2,df,'A^{-1}');
%     title(['CC =' num2str(c)]);
    
    
    
    subplot(2,3,5);
    dspect=imscale(limits.*diffsp);
    dref=imscale(limits.*CTF(nu,pixA,P).^2);
    dspect(1:nu/2+1,:)=dref(1:nu/2+1,:);
    xs=(-nu/4:nu/4-1)*df;
    imac(xs,xs,Crop(dspect,nu/2));
    xlabel('A^{-1}');
    title(['\Deltaz=' num2str(-p(1)) '   Astig=' num2str(p(2))...
        '   \theta=' num2str(180/pi*p(3))]);

    subplot(2,3,4);
    xs=(-nu/2:nu/2-1)*df;
    imac(xs,xs,dspect);
    xlabel('A^{-1}');
    title(['CC =' num2str(c)]);

    
    % Make the radial display
    subplot(2,3,6);
    rs=Radial(limits.*diffsp);
    rs=150*(rs-min(rs))/(max(rs)-min(rs));
    radCorrSpecs=[rs Radial(dref)];
    plot(freqs,radCorrSpecs);
    xlabel('A^{-1}');
    drawnow;
end;
% P.defocus=-P.defocus;  % change the polarity to match scope.
if test>1
    c
    Final_values=[P.defocus P.deltadef 180/pi*P.theta P.alpha]
end;

