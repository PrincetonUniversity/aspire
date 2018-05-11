function SetComplex
% function SetComplex
% Set the color map to allow polar representation of complex numbers.  The lowest 256
% levels are a grayscale.  Starting with element 257 we index 2^12 colors,
% based on 6 bits (64 levels) of intensity and 6 bits of color.  Thus given a complex number
% represented as (r,t), with 0<r<1 and 0<t<2*pi we have
%     index = 256+floor(r*64)+64*floor(64*t/(2*pi))
% The function imacx does this 
n=64;

r10=sin(pi/2*(1:-1/21:0));      % 22 points
% r10=1:-1/21:0;      % 22 points
r01=sin(pi/2*(0:1/21:21/22)); % 21 points
% r01=0:1/21:21/22; % 21 points
r00=zeros(1,21);
gray=0:1/255:1;     % 256 points

ared=[r10 r00 r01];
agrn=[r01 r10 r00];
ablu=[r00 r01 r10];
aint=0:1/63:1;  % 64 intensity levels

% we mix the intensities to make the brightness more constant.
mred=min(kron(aint',ared+0.1*ablu),1) ;
mgrn=min(kron(aint',0.8*agrn+0.2*ared+0.2*ablu),1);
mblu=min(kron(aint',ablu+0.2*ared),1);

vred=reshape(mred,n*n,1);
vgrn=reshape(mgrn,n*n,1);
vblu=reshape(mblu,n*n,1);

map=[gray' gray' gray'; vred vgrn vblu];
colormap(map)
