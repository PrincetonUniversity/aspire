function projs=gridMakeProjections(map, angs)
% angles is an 3 x ma matrix.  Creates na projections

n=size(map,1);
na=size(angs,2);

comp=gridMakePreComp(n,3);
P3=gridMakePaddedFT(map,'grid',comp);  % get the 3D fft in a form for slicing.

projs=single(zeros(n,n,na));  % particle stack.
for i=1:na
    P2=gridExtractPlane(P3,angs(:,i));
    projs(:,:,i)=gridRecoverRealImage(P2);
end;
