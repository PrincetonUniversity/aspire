A=rand(1000,'single');
B=rand(1000,'single');

tic;
gA=gpumatrix(A);
gB=gpumatrix(B);
for k=1:100;
    gC=gA*gB;
end
C1=gather(gC);
toc;

tic;
for k=1:100
    C2=A*B;
end
toc;


norm(C1(:)-C2(:))/norm(C2(:))