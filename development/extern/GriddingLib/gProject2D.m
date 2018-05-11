function projs=gProject2D(volume, angles)
% function planes=Project2D(volume, angles);
% Compute projections of the n x n x n volume at the Euler angles given by
% 'angles', a 3 x na array.  The output projs is an n x n x na array.  This
% function uses the GriddingLib routines.

kernelwidth=3;

[n n1 n2]=size(volume);
[n3 na]=size(angles);
if n3 ~= 3
    error('ExtractPlanes: angles argument must be 3 x na');
end;

% Make the pre-compensated 3D FFT
comp=gridMakePreComp(n,kernelwidth);
V=gridMakePaddedFT(volume,'grid',comp);

% Compute the projections
planes=zeros(n,n,na);
for i=1:na
    P=gridExtractPlane(V,angles(:,i),kernelwidth);
    projs(:,:,i)=gridRecoverRealImage(P);
end;
