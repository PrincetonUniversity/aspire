A=rand(10,'single');
B=rand(10,'single');

gA=gpumatrix(A);
gB=gpumatrix(B);
gC=gA*gB;
C1=gather(gC);

C2=A*B;

norm(C1(:)-C2(:))/norm(C2(:))