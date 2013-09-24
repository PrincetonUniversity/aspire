function p2=gridExtractPlaneE(P, angles, kernelsize)
% function p2=gridExtractPlaneE(P, angles, kernelsize);
% E version, uses standard Euler matrix convention.
% Compute the central-section plane of the ft of a volume,
% using the data structure P from the function gridMakePaddedFT.
% The plane is returned as a data structure as well, with
% p2.PadFT being np x np in size and
% in cosine-and-sine form, i.e. a real vector.
% Vectorized version.  fs 10 Feb 07

% default argument
if nargin<3
    kernelsize=3;
end;

[w1 ov]=gridMakeKaiserTable(kernelsize,'grid');
w1=w1';

% Interpolation is done here.
np=P.np;
np1=P.np1;
np1ctr=P.np1ctr;
sp1=P.sp1;
ovctr=ov/2+1;

p2=gridMakeNullFT(P.n,2);  % output structure

% Vectorized code.
plane=zeros(np*np,1);  % vector to receive the data
nw2=floor(kernelsize/2);  % half-width of kernel (=(nw-1)/2).
E=inv(EulerMatrixStandard(angles));

% Compute the central section.  The k's are implicitly set to 0.
% Source coordinates
rmax=np/2+1;
[is,js]=ndgrid(-np/2:np/2-1);
is=reshape(is,np*np,1);
js=reshape(js,np*np,1);

% Object coordinates
ip=E(1,1)*is+E(1,2)*js+np1ctr;
ip=min(sp1+np,max(ip,sp1+1));  % prevent coords from going out of bounds
ipint=round(ip);
ipfrac=floor((ip-ipint)*ov)+ovctr;

jp=E(2,1)*is+E(2,2)*js+np1ctr;
jp=min(sp1+np,max(jp,sp1+1));
jpint=round(jp);
jpfrac=floor((jp-jpint)*ov)+ovctr;

kp=E(3,1)*is+E(3,2)*js+np1ctr;
kp=min(sp1+np,max(kp,sp1+1));
kpint=round(kp);
kpfrac=floor((kp-kpint)*ov)+ovctr;

addrs=ipint+np1*((jpint-1)+np1*(kpint-1));  % 1-dim addresses in the PadFT array.
for i1=-nw2:nw2
    for j1=-nw2:nw2
        for k1=-nw2:nw2
        plane=plane+w1(ipfrac,i1+nw2+1).*w1(jpfrac,j1+nw2+1).*w1(kpfrac,k1+nw2+1)...
            .*P.PadFT(addrs+i1+(np1*(j1+np1*k1)));
        end;
    end;
end;
% Convert the vector to a 2d array
p2.PadFT(sp1+1:sp1+np,sp1+1:sp1+np)=reshape(plane,np,np);
