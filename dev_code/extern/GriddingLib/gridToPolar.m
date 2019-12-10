function p=gridToPolar(m,ntheta)
% function p=gridToPolar(m,ntheta)
% Use Gridding to convert the n x n cartesian, band-limited image m
% - or stack of images - to
% polar coordinates. The result p is n/2 x ntheta in size, where ntheta is
% by default 3 x n.


% % % test code
% n=256;
% m=zeros(n,n);
% m(:,n/4+2)=1;
% m=fuzzymask(n,2,20,2,[1,1]*(n/2+1));
% m(90:100,90:100)=1;
% m(20,20:n-20)=1;
% m=real(ifftn(fftn(m).*fftshift(fuzzymask(n,2,n/2-10,8))));
% nargin=1;

% parameters
[n ny nim]=size(m);
if nargin<2
    ntheta=3*n;
end;

mode='grid';
kernelsize=3;
gridfactor=1.25;

nw=kernelsize;  % size of interp window
nw2=floor(nw/2);  % half-width of kernel, =1 for nw=3.

np=n*gridfactor;  % padded size
np1=np+NextNiceNumber(2*nw2,7,4);  % expanded to include a border for kernel points
sp1=(np1-np)/2;  % offset into np1-sized volume.

% Make the look-up table for interpolation
[w1 ov]=gridMakeKaiserTable(kernelsize,mode);
w1=w1';

ovctr=ov/2+1;  % Center of oversampled interpolation table.

% Vectorized code.
p=zeros(n/2,ntheta,nim);

% % Source coordinates
% [is,js]=ndgrid(-np/2+sp1:np/2-1+sp1);
% is=is(:);
% js=js(:);

dt=2*pi/ntheta;

% Object coordinates: zero-based
[r t]=ndgrid(0:gridfactor:np/2-gridfactor,0:dt:2*pi-dt);
% These should be same size as p!
ip=r.*cos(t)+np1/2;  % 0-based indices
jp=r.*sin(t)+np1/2;    %
ipint=round(ip(:));
ipfrac=floor((ip(:)-ipint)*ov)+ovctr;
jpint=round(jp(:));
jpfrac=floor((jp(:)-jpint)*ov)+ovctr;

comp=gridMakePreComp(n,nw);
comp2=kron(comp,comp');

parfor i=1:nim
    m1=m(:,:,i);
    % Make the oversampled, precompensated image
    fpm=Crop(fftshift(fftn(fftshift(m1))),np);  % pad the FT of m
    pm1=Crop(fftshift(real(ifftn(fftshift(fpm.*comp2)))),np1);
    
    p1=zeros(n/2*ntheta,1);
    
    addrs=ipint+np1*(jpint)+1;  % 1-dim addresses in the PadFT array.
    for i1=-nw2:nw2
        for j1=-nw2:nw2 % Accumulate interpolated values to each point
            p1=p1+w1(ipfrac,i1+nw2+1).*w1(jpfrac,j1+nw2+1)...
                .*pm1(addrs+i1+(np1*j1));
        end;
    end;
    p1=reshape(p1,n/2,ntheta)*gridfactor^2;
    p(:,:,i)=p1;
    
end;

% % more test code
% subplot(2,2,1);
% imacs(m);
% m0=m;
% subplot(2,2,2);
% imacs(p);
% subplot(2,2,3);
% plot(sect(p));
% % pause;
% gridToRect;
