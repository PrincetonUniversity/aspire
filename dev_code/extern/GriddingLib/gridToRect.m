function m=gridToRect(p)
% function m1=gridToRect(p)
% Convert the n/2 x ntheta polar representation p into a square image m of
% dimension n x n.  We use gridding with oversampling by 1.25x.
% if p is n/2 x 1, then a circularly symmetric image is produced.
%
%  p=zeros(64,384);
%  p(1:20,:)=1;
% p(30,:)=1;

[n2 ntheta]=size(p);
n=2*n2;

if ntheta<2
    ntheta=6*n2;
    p=repmat(p,1,ntheta);
end;

kernelsize=3;
gridfactor=1.25;

nw=kernelsize;  % size of interp window
nw2=floor(nw/2);  % half-width of kernel, =1 for nw=3.

np=n*gridfactor;  % padded size
np1=np+NextNiceNumber(2*nw2,7,4);  % expanded to include a border for kernel points
sp1=(np1-np)/2;  % offset into np1-sized volume.
[w1 ov]=gridMakeKaiserTable(kernelsize,'grid');
% w1(:,1025)=w1(:,1024);  % extend the table by 1.
% Get the extra-padded plane.
mp1=zeros(np1,np1);
cp1=np1/2+1;  % Center of the padded plane.

ovctr=ov/2+1;
square=zeros(nw,nw);
% loop over r and t, the coordinates of the input p.
dt=2*pi/ntheta;
for r=1:n2
    if r>1
        wt=dt/2*((r-.5)^2-(r-1.5)^2); % pi/nthteta*r^2 of annulus from r-.5 to r+.5
    else
        wt=dt/8/sqrt(2);  % dt/2*(.5)^2.
        % I have no idea what the sqrt(2) is doing....d
    end;
    
    for t=1:ntheta
        r0=(r-1)*gridfactor;
        ip=r0*cos((t-1)*dt)+cp1;
        jp=r0*sin((t-1)*dt)+cp1;
        ipint=round(ip);
        ipfrac=floor((ip-ipint)*ov)+ovctr;
        ip0=ipint;
        jpint=round(jp);
        jpfrac=floor((jp-jpint)*ov)+ovctr;
        jp0=jpint;
        % we compute the interpolated function on a square surrounding (ip jp)
        % one input point at a time.
        square=p(r,t)*w1(:,ipfrac)*w1(:,jpfrac)';
        % Add the square into the plane.
        mp1(ip0-nw2:ip0+nw2,jp0-nw2:jp0+nw2)=...
            mp1(ip0-nw2:ip0+nw2,jp0-nw2:jp0+nw2)+wt*square;
    end; % t
end; %r

% transform and correct
fmp=fftn(Crop(mp1,np));  % strip out-of-bounds points and compute ft
comp=gridMakePreComp(n,nw);
%   comp=1;
% Multiply by the compensation function, ift, center and crop to n x n.
fm=Crop(fftshift(fmp).*kron(comp,comp'),n).*fuzzymask(n,2,n/2-2,4);
m=real(ifftn(fftshift(fm)));
% 
% subplot(2,2,3); imacs(p);
% subplot(2,2,4); imacs(m);
% subplot(2,2,2);
% plot([sect(m) sect(m0)]);
% plot([sect(m0) sect(m)]);