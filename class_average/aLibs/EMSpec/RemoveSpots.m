function [spq qmask]=RemoveSpots(sp,minr,ThreshSD,display)
% function [spc pmask]=RemoveSpots(sp,minr)
% Simple program to remove crystal periodicity from a centered power spectrum
% We ignore the spectrum out to a radius of minr.  This can be computed as
% f*(n*pixA)/58.2 where f is a factor (e.g. 0.8) to allow for error.  Here
% 58.2A is the lattics constant of the streptavidin crystal and pixA is
% Angstroms per pixel.
minr=ceil(minr);
% n=size(fd2,1);
% sp=abs(fd2).^2;
% minr=100;


sz0=size(sp);
n0=max(sz0);

if nargin<3
    ThreshSD=6;
end;
if nargin<4
    display=0;
end;
ThreshSDs=[1 1 1]*ThreshSD;
% end;

% spq=Downsample(sp,n0);  % working spectrum, made square if necessary
spq=sp;

for i=1:numel(ThreshSDs)
    psp=ToPolar(spq);
    [nr nt]=size(psp);
    psp(1:minr,:)=0;
    
    mes=median(psp')';  % std at each circle.
    sds=std(psp')';
    mex=repmat(mes,1,nt);
    sdx=repmat(sds,1,nt);
    sdx(1:minr,:)=max(sp(:));
    mec=ToRect(mex);
    sdc=ToRect(sdx);
    nsp=(sp-mec)./sdc;
    excl=single(nsp>ThreshSDs(i));
    
    qmask=GaussFilt(excl,.1)<.1;
    spq=spq.*qmask;
    
    if display
        subplot(2,3,i);
        imacs(qmask);
        % imacs(sdc);
        drawnow;
    end;
    % psp2=psp;
end;

% pmask=Downsample(qmask,sz0); % change the aspect ratio of the mask if necessary.
pmask=qmask;
spc=sp.*pmask;

if display
    subplot(2,3,i+1);
    imacs((spq).^.1);
    
    subplot(2,3,i+2);
    semilogy(Radial(spq));
end;

