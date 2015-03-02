function OUTPUT=fastrotate3d_gpu(INPUT,R,precision)
%FASTROTATE3D_GPU Rotate a 3D volume by a given rotation matrix.
%
% GPU implementation of fastrotate3d.
% See fastrotate3d for information.
%
%Yoel Shkolnisky, April 2014.

if ~exist('precision','var')
    precision='single';
end


[psi,theta,phi]=rot2xyz(R);

psid=psi*180/pi;
thetad=theta*180/pi;
phid=phi*180/pi;

[SzX, SzY, SzZ] =size(INPUT);

if strcmpi('precision','single')
    INPUT=single(INPUT);
end

gINPUT=gpuArray(INPUT); % Set precision
gOUTPUT=gpuArray.zeros(size(INPUT),precision);


%tmp=fastrotate3x(INPUT,psid);
M=fastrotateprecomp(SzY,SzZ,psid);

if strcmpi('precision','single')
    M.Mx=single(M.Mx);
    M.My=single(M.My);
end
gMx=gpuArray(M.Mx);
gMy=gpuArray(M.My);


for k=1:SzX
    gim=squeeze(gINPUT(:,k,:));
    gOUTPUT(:,k,:)=fastrotate3daux_gpu(gim,gMx,gMy,M.mult90,precision);    
end


%tmp=fastrotate3y(tmp,thetad);
M=fastrotateprecomp(SzX,SzZ,-thetad);

if strcmpi('precision','single')
    M.Mx=single(M.Mx);
    M.My=single(M.My);
end
gMx=gpuArray(M.Mx);
gMy=gpuArray(M.My);

for k=1:SzY
    gim=squeeze(gOUTPUT(k,:,:));
    gOUTPUT(k,:,:)=fastrotate3daux_gpu(gim,gMx,gMy,M.mult90,precision);       
end

%OUTPUT=fastrotate3z(tmp,phid);

M=fastrotateprecomp(SzX,SzY,-phid);

if strcmpi('precision','single')
    M.Mx=single(M.Mx);
    M.My=single(M.My);
end
gMx=gpuArray(M.Mx);
gMy=gpuArray(M.My);

for k=1:SzZ
    gim=squeeze(gOUTPUT(:,:,k));
    gOUTPUT(:,:,k)=fastrotate3daux_gpu(gim,gMx,gMy,M.mult90,precision);    
end

OUTPUT=double(gather(gOUTPUT));