function P1=gridInsertPlanePar(p2, P, angles, kernelsize, mode)
% Experimental parallelized version of gridInsertPlane--doesn't work
% because of complex indexing at line 73.
% function P1=gridInsertPlane(p2, P, angles[, kernelsize, mode]);
% Add the central-section plane to the ft of a volume,
% using the 2d data structure p2 and the 3D structure P
% as obtained from the function gridMakePaddedFT.
% The plane p2.PadFT is assumed to be
% in cosine-and-sine form, i.e. a real 2d array.
% angles are the Euler angles.
% mode is either 'sinc' or 'grid'; default is 'grid'
% default kernelsize is 3.
%
% Vectorized code doesn't work, this has the "scalar" code.

% default arguments
if nargin<5
    mode='grid';
end;
if nargin<4
    kernelsize=3;
end;

nw=kernelsize;  % size of interp window


[w1 ov]=gridMakeKaiserTable(kernelsize,mode);


% Interpolation is done here.
%"scalar" code.

E=inv(EulerMatrix(angles));  % use the inverse transformation.
P1=P;  % Copy the original volume.
np=P.np;
np1=P.np1;
np1ctr=P.np1ctr;
sp1=P.sp1;

if p2.np1 ~= np1
    error(['gridInsertPlane: plane has wrong dimension, ', num2str(P.np1)]);
end;

% Get the extra-padded plane.
pplane=p2.PadFT;
r=np/2+1;  % radius of plane region to sample
% Do the interpolated insertion
centerp=ones(3,1)*np1ctr;
centerp2=centerp(1:2);
nw2=floor(nw/2);  % half-width of kernel (=(nw-1)/2).
ovctr=ov/2+1;
% loop over i,j,k the coordinates of the plane's system
pft=P1.PadFT;
parfor i=np1ctr-np/2:np1ctr+np/2-1
    jspan=floor(sqrt(r^2-(i-np1ctr)^2));
    for j=np1ctr-jspan:np1ctr+jspan
        k=np1ctr;
        p=E*([i;j;k]-centerp)+centerp;
        pint=round(p);
        pfrac=floor((p-pint)*ov)+ovctr;
         % new plane coordinates
        i1=pint(1);
        j1=pint(2);
        k1=pint(3);
         % we compute the interpolated function on a cube surrounding (i1 j1 k1)
         % one plane at a time.
        insplane=pplane(i,j)*w1(:,pfrac(1))*w1(:,pfrac(2))';
cube=zeros(nw,nw,nw);
        for k2=1:nw
            cube(:,:,k2)=insplane*w1(k2,pfrac(3));
        end;
         % Add the cube into the volume.
        pft(i1-nw2:i1+nw2,j1-nw2:j1+nw2,k1-nw2:k1+nw2)...
            =pft(i1-nw2:i1+nw2,j1-nw2:j1+nw2,k1-nw2:k1+nw2)+cube;
    end; % j
end; %i
P1.PadFT=pft;